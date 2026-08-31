#!/usr/bin/env bash
# ============================================================================
# 胚系（Germline）GPU 主流程 第一步
#   工具：NVIDIA Parabricks  pbrun deepvariant_germline
#   一条命令完成：FASTQ -> fq2bam(比对+排序+markdup+BQSR) -> DeepVariant -> GVCF
#
#   注意（与 CPU 标准流程的差异，务必知悉）：
#   - 本流程使用 DeepVariant（深度学习 caller），不使用 GATK HaplotypeCaller。
#   - DeepVariant 路线使用模型输出的变异评分，不沿用 HaplotypeCaller 路线的 GATK VQSR。
#     本脚本在联合分型后按 GQ/DP 做基因型级过滤；阈值需根据队列与验证结果调整。
#   - 参考基因组与 BQSR known-sites 均与 CPU 标准流程对齐，两套流程唯一的本质区别是
#     caller（HaplotypeCaller vs DeepVariant）与联合分型工具（GenotypeGVCFs vs GLnexus）。
#   - 本步骤输出每个样本的 GVCF（--gvcf），供第二步 GLnexus 联合基因分型使用。
#     GLnexus 用于 DeepVariant GVCF 的 cohort-level joint genotyping；
#     CPU 路线则使用 GenomicsDBImport + GenotypeGVCFs。
#
#   投递方式（GPU host without a scheduler，用 nohup 后台串行跑）：
#     cd /path/to/germline-project
#     nohup bash scripts/gpu/01_deepvariant.sh /path/to/germline-project/samples.list \
#       > germline_gpu_step1.main.log 2>&1 &
#
#   如需改路径/参数，可用环境变量覆盖，例如：
#     DATASET_DIR=/path/to/germline-project/data \
#     BASE_OUT=/path/to/germline-project \
#     REF_ROOT=/path/to/reference \
#     SEQ_MODE=wgs \
#     nohup bash scripts/gpu/01_deepvariant.sh samples.list \
#       > step1.main.log 2>&1 &
# ============================================================================
set -euo pipefail

###########################
# 加载环境变量（取 fastp 等，如不需要可忽略缺失）
###########################
ENV_FILE="${ENV_FILE:-}"
if [[ -s "$ENV_FILE" ]]; then
  source "$ENV_FILE"
  echo "✓ 已加载环境文件: $ENV_FILE"
else
  echo "警告：环境文件不存在: $ENV_FILE（若不依赖其中工具可忽略）" >&2
fi

###########################
# 配置参数（允许环境变量覆盖）
###########################
# 输出根目录
BASE_OUT="${BASE_OUT:-${PWD}/work}"

# 样本列表（第一个位置参数；每行一个样本名）
SAMPLE_LIST="${1:-${BASE_OUT}/samples.list}"

# 数据与参考路径（宿主机路径）
DATASET_DIR="${DATASET_DIR:-${BASE_OUT}/data}"
REF_ROOT="${REF_ROOT:-${BASE_OUT}/reference}"

# Parabricks 容器镜像示例；请在正式运行前按目标环境确认并固定版本
IMAGE="${IMAGE:-nvcr.io/nvidia/clara/clara-parabricks:4.5.1-1}"

# 测序类型：wgs / wes（影响 fastp 数据量判断与后续 GLnexus 配置；deepvariant 命令本身一致）
SEQ_MODE="${SEQ_MODE:-wgs}"

# Read group 的 LB 字段
RG_LB="${RG_LB:-WGS}"

# 是否在 deepvariant 前先用 fastp 单独质控（1=做，0=跳过直接喂原始 FASTQ）
# 说明：fq2bam 内部不做接头/质量修剪，建议保留 fastp 质控这一步
DO_FASTP="${DO_FASTP:-1}"
FASTP="${RESEQ_PREFIX:-/usr/local}/bin/fastp"

# 容器内路径（固定挂载点，勿改）
REF="/refdir/GRCh38.fasta"
# BQSR known-sites —— 与 CPU 路线使用同一组示例资源
# 正式运行前请根据参考基因组版本、validated reference bundle 和本地验证结果确认资源。
KNOWN1="/refdir/dbsnp_146.hg38.vcf.gz"
KNOWN2="/refdir/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
KNOWN3="/refdir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
KNOWN4="/refdir/1000G_omni2.5.hg38.vcf.gz"

###########################
# 准备
###########################
sed -i 's/\r$//' "${SAMPLE_LIST}" 2>/dev/null || true
[[ -s "${SAMPLE_LIST}" ]] || { echo "[FATAL] 找不到/为空的样本列表: ${SAMPLE_LIST}" >&2; exit 1; }
[[ -s "${REF_ROOT}/GRCh38.fasta" ]] || { echo "[FATAL] 找不到参考: ${REF_ROOT}/GRCh38.fasta" >&2; exit 1; }

mkdir -p "${BASE_OUT}/01_deepvariant/log" \
         "${BASE_OUT}"/{BQSR,gvcf,vcf,clean_fq,qc_reports,tmp}

###########################
# 耗时/状态记录辅助函数
#   每个样本写自己的 timing.log（${BASE_OUT}/01_deepvariant/log/<sample>.timing.log）
#   便于批量统计进度，例如：
#     grep -l "END.*deepvariant" 01_deepvariant/log/*.timing.log | wc -l   # 多少样本跑完
#     grep -L "ALL DONE"        01_deepvariant/log/*.timing.log            # 哪些样本没跑完
###########################
fmt_time() { local s=$1; printf '%dh%02dm%02ds' $((s/3600)) $(((s%3600)/60)) $((s%60)); }

TIMING=""          # 当前样本的 timing.log 路径（进入样本循环时设置）
STEP_START_TS=0

step_start() {  # $1=步骤名
  STEP_START_TS=$(date +%s)
  echo "[$(date '+%F %T')] 开始 ${SAMPLE_CUR} $1..."
  printf '[%s] START  %s\n' "$(date '+%F %T')" "$1" >> "${TIMING}"
}
step_end() {    # $1=步骤名
  local end_ts elapsed; end_ts=$(date +%s); elapsed=$((end_ts - STEP_START_TS))
  echo "[$(date '+%F %T')] 完成 ${SAMPLE_CUR} $1，耗时 $(fmt_time ${elapsed})"
  printf '[%s] END    %-14s elapsed=%-7s (%s)\n' \
    "$(date '+%F %T')" "$1" "${elapsed}s" "$(fmt_time ${elapsed})" >> "${TIMING}"
}

echo "=========================================="
echo "胚系 GPU 第一步：pbrun deepvariant_germline (FASTQ -> GVCF)"
echo "样本列表 : ${SAMPLE_LIST}"
echo "镜像     : ${IMAGE}"
echo "参考     : ${REF_ROOT}/GRCh38.fasta"
echo "测序类型 : ${SEQ_MODE}    fastp质控: ${DO_FASTP}"
echo "=========================================="

###########################
# 逐样本处理（串行；GPU host without a scheduler）
###########################
while read -r sample; do
  [[ -z "${sample}" ]] && continue

  SAMPLE_CUR="${sample}"
  echo ""
  echo ">>> 处理样本: ${sample}    开始: $(date)"

  raw_fq1="${DATASET_DIR}/${sample}/${sample}_f1.fastq.gz"
  raw_fq2="${DATASET_DIR}/${sample}/${sample}_r2.fastq.gz"
  [[ -f "${raw_fq1}" ]] || { echo "[WARN] 缺少 ${raw_fq1}，跳过 ${sample}" >&2; continue; }
  [[ -f "${raw_fq2}" ]] || { echo "[WARN] 缺少 ${raw_fq2}，跳过 ${sample}" >&2; continue; }

  out_bam="${BASE_OUT}/BQSR/${sample}.BQSR.bam"
  out_gvcf="${BASE_OUT}/gvcf/${sample}.g.vcf.gz"
  recal="${BASE_OUT}/BQSR/${sample}.BQSR_REPORT.txt"
  log_file="${BASE_OUT}/01_deepvariant/log/${sample}.log"
  TIMING="${BASE_OUT}/01_deepvariant/log/${sample}.timing.log"
  : > "${TIMING}"   # 每次运行该样本重置其 timing（断点续跑下方会先判断 GVCF 是否已存在）

  # 断点续跑：GVCF 已存在则跳过
  if [[ -s "${out_gvcf}" ]]; then
    echo ">>> 跳过（GVCF 已存在）: ${sample}"
    printf '[%s] SKIP   GVCF already exists\n' "$(date '+%F %T')" >> "${TIMING}"
    continue
  fi

  #############################################
  # 步骤1：fastp 质控（可选）
  #############################################
  if [[ "${DO_FASTP}" == "1" ]]; then
    clean_fq1="${BASE_OUT}/clean_fq/${sample}_clean_f1.fq.gz"
    clean_fq2="${BASE_OUT}/clean_fq/${sample}_clean_r2.fq.gz"
    if [[ -s "${clean_fq1}" && -s "${clean_fq2}" ]]; then
      echo ">>> [1/2] fastp（已完成，跳过）"
    else
      step_start fastp
      "${FASTP}" -l 50 -w 16 --compression 2 -5 -3 \
        -i "${raw_fq1}" -I "${raw_fq2}" \
        -o "${clean_fq1}" -O "${clean_fq2}" \
        -j "${BASE_OUT}/qc_reports/${sample}.fastp.json" \
        -h "${BASE_OUT}/qc_reports/${sample}.fastp.html" \
        >> "${log_file}" 2>&1
      [[ -s "${clean_fq1}" && -s "${clean_fq2}" ]] || { echo "[WARN] fastp 失败，跳过 ${sample}" >&2; continue; }
      step_end fastp
    fi
    in_fq1="/workdir/clean_fq/${sample}_clean_f1.fq.gz"
    in_fq2="/workdir/clean_fq/${sample}_clean_r2.fq.gz"
  else
    # 不做 fastp，直接把原始 FASTQ 软链到 clean_fq 以统一容器内路径
    in_fq1="/workdir/clean_fq/${sample}_clean_f1.fq.gz"
    in_fq2="/workdir/clean_fq/${sample}_clean_r2.fq.gz"
    ln -sf "${raw_fq1}" "${BASE_OUT}/clean_fq/${sample}_clean_f1.fq.gz"
    ln -sf "${raw_fq2}" "${BASE_OUT}/clean_fq/${sample}_clean_r2.fq.gz"
  fi

  #############################################
  # 步骤2：Parabricks deepvariant_germline
  #   一条命令：比对+BQSR(fq2bam) + DeepVariant 出 GVCF
  #############################################
  echo ">>> [2/2] pbrun deepvariant_germline (GPU)..."
  step_start deepvariant
  docker run --gpus all --rm \
    -u "$(id -u)":"$(id -g)" \
    --workdir /workdir \
    -v "${REF_ROOT}":/refdir \
    -v "${BASE_OUT}":/workdir \
    "${IMAGE}" pbrun deepvariant_germline \
      --ref "${REF}" \
      --in-fq "${in_fq1}" "${in_fq2}" \
              "@RG\tID:${sample}\tLB:${RG_LB}\tPL:Illumina\tSM:${sample}\tPU:${sample}" \
      --knownSites "${KNOWN1}" \
      --knownSites "${KNOWN2}" \
      --knownSites "${KNOWN3}" \
      --knownSites "${KNOWN4}" \
      --out-bam "/workdir/BQSR/${sample}.BQSR.bam" \
      --out-recal-file "/workdir/BQSR/${sample}.BQSR_REPORT.txt" \
      --out-variants "/workdir/gvcf/${sample}.g.vcf.gz" \
      --gvcf \
      --tmp-dir "/workdir/tmp" \
      --num-streams-per-gpu 4 \
      --gpusort \
      --gpuwrite \
      >> "${log_file}" 2>&1
  step_end deepvariant

  if [[ -s "${out_gvcf}" && -s "${out_bam}" ]]; then
    echo ">>> 完成: ${sample}"
    echo "    GVCF: ${out_gvcf}"
    echo "    BAM : ${out_bam}"
    printf '[%s] ALL DONE  %s\n' "$(date '+%F %T')" "${sample}" >> "${TIMING}"
    # 清理质控临时 FASTQ（仅当 GVCF 与 BAM 都正常产出）
    if [[ "${DO_FASTP}" == "1" ]]; then
      rm -f "${BASE_OUT}/clean_fq/${sample}_clean_f1.fq.gz" \
            "${BASE_OUT}/clean_fq/${sample}_clean_r2.fq.gz"
    fi
  else
    echo "[WARN] 失败或输出不完整: ${sample}（详见 ${log_file}）" >&2
    printf '[%s] FAILED    %s（输出不完整，详见 %s）\n' "$(date '+%F %T')" "${sample}" "${log_file}" >> "${TIMING}"
  fi

  echo ">>> 结束: $(date)"
  echo "---"
done < "${SAMPLE_LIST}"

echo ""
echo "=========================================="
echo "胚系 GPU 第一步完成！"
echo "  GVCF 目录   : ${BASE_OUT}/gvcf"
echo "  BAM  目录   : ${BASE_OUT}/BQSR"
echo "  每样本日志  : ${BASE_OUT}/01_deepvariant/log/<sample>.log"
echo "  每样本耗时  : ${BASE_OUT}/01_deepvariant/log/<sample>.timing.log"
echo ""
echo "  进度排查示例："
echo "    跑完的样本数 : grep -l 'ALL DONE' ${BASE_OUT}/01_deepvariant/log/*.timing.log | wc -l"
echo "    没跑完的样本 : grep -L 'ALL DONE' ${BASE_OUT}/01_deepvariant/log/*.timing.log"
echo "    各样本耗时   : grep 'END' ${BASE_OUT}/01_deepvariant/log/*.timing.log"
echo ""
echo "  下一步：运行 scripts/gpu/02_glnexus.sh 做联合基因分型"
echo "=========================================="
