#!/usr/bin/env bash
# ========================================
# 胚系流程 第一步 - fastp + BWA-MEM + markdup + BQSR
# ========================================
# 本步骤为“按样本”投递：每个样本生成一个子脚本，独立运行 fastp → BWA → 排序 → 标记重复 → BQSR
# 输出 BQSR 后的 BAM，作为 step2（HaplotypeCaller）的输入
#
# 通过传参，投递任务命令，可参考如下：
#
# DATASET_DIR=/path/to/germline-project/data \
# BASE_OUT=/path/to/germline-project \
# REF_ROOT=/path/to/reference \
# ENV_FILE=/path/to/env.sh \
# GATK=gatk \
# CORE=16 \
# MEM=96 \
# RG_LB=WES \
# SUBMIT=0 \
# bash scripts/cpu/01_preprocess_bqsr.sh /path/to/germline-project/samples.list

set -euo pipefail

###########################
# 配置参数（允许 env 覆盖）
###########################
# 输出根目录
BASE_OUT="${BASE_OUT:-${PWD}/work}"

DATASET_DIR="${DATASET_DIR:-${BASE_OUT}/data}"
REF_ROOT="${REF_ROOT:-${BASE_OUT}/reference}"

QUEUE="${QUEUE:-compute}"

# 资源配置
CORE="${CORE:-${STEP1_CORE:-8}}"
MEM="${MEM:-${STEP1_MEM:-64}}"
STEP1_CORE="${STEP1_CORE:-$CORE}"
STEP1_MEM="${STEP1_MEM:-$MEM}"
WALLTIME="${WALLTIME:-1000:00:00}"

# 是否提交任务：1=生成子脚本后用 sbatch 提交；0=只生成子脚本，不提交
# 在没有 Slurm/sbatch 的服务器上，请用 SUBMIT=0
SUBMIT="${SUBMIT:-1}"

# 环境文件 / 工具
ENV_FILE="${ENV_FILE:-}"
GATK="${GATK:-gatk}"

# 参考文件（也允许外部覆盖）
REF_FASTA="${REF_FASTA:-${REF_ROOT}/GRCh38.fasta}"
# BQSR known-sites（示例使用4个常见 known-sites 资源，可按环境调整）
KNOWN1="${KNOWN1:-${REF_ROOT}/dbsnp_146.hg38.vcf.gz}"
KNOWN2="${KNOWN2:-${REF_ROOT}/1000G_phase1.snps.high_confidence.hg38.vcf.gz}"
KNOWN3="${KNOWN3:-${REF_ROOT}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz}"
KNOWN4="${KNOWN4:-${REF_ROOT}/1000G_omni2.5.hg38.vcf.gz}"

# Read group 的 LB 字段，WES/WGS 等
RG_LB="${RG_LB:-WES}"

# --- 覆盖样本清单 / 单样本直投 ---
# 如果传进来的是一个存在的文件，就用这个文件当 SAMPLE_LIST；
# 如果传进来的是一个/多个样本名，就临时写成一个清单文件；
# 如果什么都不传，就用 DEFAULT_LIST 当 SAMPLE_LIST。
DEFAULT_LIST="${DEFAULT_LIST:-${BASE_OUT}/samples.list}"

if [[ $# -gt 0 ]]; then
  if [[ -f "$1" ]]; then
    SAMPLE_LIST="$1"
  else
    SAMPLE_LIST="$(mktemp)"
    printf '%s\n' "$@" > "$SAMPLE_LIST"
    trap 'rm -f "$SAMPLE_LIST"' EXIT
  fi
else
  SAMPLE_LIST="$DEFAULT_LIST"
fi

sed -i 's/\r$//' "$SAMPLE_LIST"
[[ -s "$SAMPLE_LIST" ]] || { echo "[FATAL] 样本清单为空：$SAMPLE_LIST" >&2; exit 1; }

mkdir -p "${BASE_OUT}/01_bwa_bqsr/sh" \
         "${BASE_OUT}/tmp" \
         "${BASE_OUT}/raw_bam" \
         "${BASE_OUT}/marked" \
         "${BASE_OUT}/BQSR"

#############################################
# 第一步：fastp + BWA + markdup + BQSR
#############################################
echo "========================================"
echo ">> 胚系 第一步：fastp + BWA-MEM + markdup + BQSR"
echo "========================================"

while read -r sample; do
  [[ -z "${sample}" ]] && continue

  shp="${BASE_OUT}/01_bwa_bqsr/sh/${sample}.sh"
  cat > "${shp}" <<'EOSH'
#!/usr/bin/env bash
#SBATCH -p __QUEUE__
#SBATCH -J g_step1___SAMPLE__
#SBATCH -c __CORE__
#SBATCH --mem=__MEM__G
#SBATCH -t __WALLTIME__
#SBATCH -o __BASE__/01_bwa_bqsr/sh/__SAMPLE__.out
#SBATCH -e __BASE__/01_bwa_bqsr/sh/__SAMPLE__.err


set -euo pipefail

[[ -n "__ENV_FILE__" && -s "__ENV_FILE__" ]] && source "__ENV_FILE__"

FASTP="${RESEQ_PREFIX}/bin/fastp"
BWA="${RESEQ_PREFIX}/bin/bwa"
SAMTOOLS="${RESEQ_PREFIX}/bin/samtools"
SAMBAMBA="${RESEQ_PREFIX}/bin/sambamba"
GATK="__GATK__"

index="__REF_FASTA__"

# 输入文件
f1="__DATASET_DIR__/__SAMPLE__/__SAMPLE___f1.fastq.gz"
f2="__DATASET_DIR__/__SAMPLE__/__SAMPLE___r2.fastq.gz"

[[ ! -f "${f1}" ]] && { echo "错误：找不到 ${f1}" >&2; exit 1; }
[[ ! -f "${f2}" ]] && { echo "错误：找不到 ${f2}" >&2; exit 1; }

# 临时文件
fq1="__BASE__/tmp/__SAMPLE___clean_1.fq.gz"
fq2="__BASE__/tmp/__SAMPLE___clean_2.fq.gz"
json="__BASE__/raw_bam/__SAMPLE__.fastp.json"
html="__BASE__/raw_bam/__SAMPLE__.fastp.html"
raw_bam="__BASE__/raw_bam/__SAMPLE__.bam"
marked_bam="__BASE__/marked/__SAMPLE__.marked.bam"
final_bam="__BASE__/BQSR/__SAMPLE__.BQSR.bam"

echo "[$(date)] 开始 fastp 质控..."
"${FASTP}" -l 50 -w __CORE__ -5 -3 \
  -i "${f1}" -I "${f2}" \
  -o "${fq1}" -O "${fq2}" \
  -j "${json}" \
  -h "${html}"

echo "[$(date)] 开始 BWA 比对 + 排序..."
"${BWA}" mem -t __CORE__ -T 0 \
  -R "@RG\tID:__SAMPLE__\tSM:__SAMPLE__\tLB:__RG_LB__\tPL:Illumina\tPU:__SAMPLE__" \
  "${index}" "${fq1}" "${fq2}" \
  | "${SAMTOOLS}" sort -@ __CORE__ -o "${raw_bam}"

echo "[$(date)] 开始 MarkDuplicates..."
# 去重复这一步极其容易出现文件句柄不足问题

# 1) 提升可打开文件数
cur=$(ulimit -n  2>/dev/null || echo 1024)
hard=$(ulimit -Hn 2>/dev/null || echo 65535)
[[ "$hard" =~ ^[0-9]+$ ]] || hard=65535
target=$(( hard < 65535 ? hard : 65535 ))
if [ "$cur" -lt "$target" ]; then
  ulimit -n "$target" 2>/dev/null || echo "[WARN] 无法提升 ulimit -n，当前=$cur 硬上限=$hard"
fi
echo "[INFO] ulimit -n = $(ulimit -n)"

# 2) 线程与临时目录
threads=__CORE__
tmpdir="${TMPDIR:-__BASE__/tmp}"
mkdir -p "$tmpdir"

# 3) 去重复（失败降线程重试一次，仍失败则退出）
if ! "${SAMBAMBA}" markdup \
      -t "$threads" \
      --tmpdir="$tmpdir" \
      "${raw_bam}" "${marked_bam}"; then
  echo "[WARN] markdup 失败，尝试用更少线程重试 (-t 4)…"
  "${SAMBAMBA}" markdup \
      -t 4 \
      --tmpdir="$tmpdir" \
      "${raw_bam}" "${marked_bam}" \
  || { echo "[FATAL] markdup 连续失败，样本=__SAMPLE__"; exit 1; }
fi

echo "[$(date)] 为 marked BAM 建立索引..."
"${SAMTOOLS}" index -@ __CORE__ "${marked_bam}"

echo "[$(date)] 开始 BaseRecalibrator..."
"${GATK}" --java-options "-Xms__MEM__G -Xmx__MEM__G -XX:ParallelGCThreads=__CORE__" \
  BaseRecalibrator \
  -R "${index}" \
  --known-sites "__KNOWN1__" \
  --known-sites "__KNOWN2__" \
  --known-sites "__KNOWN3__" \
  --known-sites "__KNOWN4__" \
  -I "${marked_bam}" \
  -O "__BASE__/BQSR/__SAMPLE___recal_data.table1"

echo "[$(date)] 开始 ApplyBQSR..."
"${GATK}" --java-options "-Xms__MEM__G -Xmx__MEM__G -XX:ParallelGCThreads=__CORE__" \
  ApplyBQSR \
  -R "${index}" \
  -I "${marked_bam}" \
  --bqsr-recal-file "__BASE__/BQSR/__SAMPLE___recal_data.table1" \
  -O "${final_bam}"

echo "[$(date)] 索引最终 BAM..."
"${SAMTOOLS}" index -@ __CORE__ "${final_bam}"

# 仅当最终 BAM 及其索引都存在时再删除中间 FASTQ
if [[ -s "${final_bam}" && -s "${final_bam}.bai" ]]; then
  skip_clean=
  case "${fq1}" in
    __BASE__/*) : ;;
    *) echo "安全检查：fq1 不在输出目录，跳过删除"; skip_clean=1 ;;
  esac
  case "${fq2}" in
    __BASE__/*) : ;;
    *) echo "安全检查：fq2 不在输出目录，跳过删除"; skip_clean=1 ;;
  esac

  if [[ -z "${skip_clean:-}" ]]; then
    echo "[$(date)] 清理临时 FASTQ: ${fq1} ${fq2}"
    rm -f "${fq1}" "${fq2}"
  else
    echo "[$(date)] 跳过清理 FASTQ（安全检查未通过）"
  fi
fi

echo "[$(date)] 第一步完成：__SAMPLE__"
EOSH

  sed -i \
    -e "s#__QUEUE__#${QUEUE}#g" \
    -e "s#__SAMPLE__#${sample}#g" \
    -e "s#__RG_LB__#${RG_LB}#g" \
    -e "s#__CORE__#${CORE}#g" \
    -e "s#__MEM__#${MEM}#g" \
    -e "s#__WALLTIME__#${WALLTIME}#g" \
    -e "s#__BASE__#${BASE_OUT}#g" \
    -e "s#__ENV_FILE__#${ENV_FILE}#g" \
    -e "s#__REF_FASTA__#${REF_FASTA}#g" \
    -e "s#__DATASET_DIR__#${DATASET_DIR}#g" \
    -e "s#__KNOWN1__#${KNOWN1}#g" \
    -e "s#__KNOWN2__#${KNOWN2}#g" \
    -e "s#__KNOWN3__#${KNOWN3}#g" \
    -e "s#__KNOWN4__#${KNOWN4}#g" \
    -e "s#__GATK__#${GATK:-gatk}#g" \
    "${shp}"

  chmod +x "${shp}"

  if [[ "${SUBMIT}" == "1" ]]; then
    if command -v sbatch >/dev/null 2>&1; then
      jid=$(sbatch "${shp}" | awk '{print $4}')
      echo "   提交: ${sample} -> Job ID: ${jid}"
    else
      echo "[FATAL] SUBMIT=1 但当前环境找不到 sbatch。" >&2
      echo "        the current environment请用 SUBMIT=0 只生成子脚本，然后用 nohup bash 子脚本运行。" >&2
      exit 2
    fi
  else
    echo "   只生成不提交: ${sample} -> ${shp}"
  fi
done < "${SAMPLE_LIST}"

echo ""
if [[ "${SUBMIT}" == "1" ]]; then
  echo ">> 第一步所有任务已提交完成！"
  echo ">> 等待所有样本的 BWA+BQSR 完成后，再运行 step2（HaplotypeCaller）"
else
  echo ">> 第一步子脚本已生成完成，但没有提交。"
  echo ">> 你可以用: nohup bash <子脚本路径> > <日志路径> 2>&1 &"
fi
echo "========================================"
