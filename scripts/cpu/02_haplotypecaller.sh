#!/usr/bin/env bash
# ========================================
# 胚系流程 第二步 - HaplotypeCaller (GVCF 模式)
# ========================================
# 本步骤为“按样本”投递：每个样本对 BQSR 后的 BAM 单独运行 HaplotypeCaller，
# 使用 -ERC GVCF 输出包含所有位点信息的 GVCF 文件，作为 step3 联合基因分型的输入。
#
# 投递示例：
# BASE_OUT=/path/to/germline-project \
# REF_ROOT=/path/to/reference \
# ENV_FILE=/path/to/env.sh \
# GATK=gatk \
# CORE=8 \
# MEM=32 \
# SUBMIT=0 \
# bash scripts/cpu/02_haplotypecaller.sh /path/to/germline-project/samples.list

set -euo pipefail
[ -n "${BASH_VERSION:-}" ] || { echo "请用: bash $0"; exit 2; }

###########################
# 配置参数（允许 env 覆盖）- 对齐 step1 习惯
###########################
# 输出根目录（必须与 step1 的 BASE_OUT 保持一致，否则找不到 BQSR bam）
BASE_OUT="${BASE_OUT:-${PWD}/work}"

# 参考与环境
ENV_FILE="${ENV_FILE:-}"
REF_ROOT="${REF_ROOT:-${BASE_OUT}/reference}"
GATK="${GATK:-gatk}"

# 资源
QUEUE="${QUEUE:-compute}"
CORE="${CORE:-${STEP2_CORE:-8}}"
MEM="${MEM:-${STEP2_MEM:-32}}"
WALLTIME="${WALLTIME:-1000:00:00}"

# Java 内存：与申请内存保持一致
JAVA_MEM_G="${JAVA_MEM_G:-${MEM}}"

# 参考文件（允许外部覆盖）
REF_FASTA="${REF_FASTA:-${REF_ROOT}/GRCh38.fasta}"
DBSNP="${DBSNP:-${REF_ROOT}/dbsnp_146.hg38.vcf.gz}"
# 分析区间：WGS 用 wgs_calling_regions；WES 请用目标区域 interval_list 覆盖
INTERVALS="${INTERVALS:-${REF_ROOT}/hg38_v0_wgs_calling_regions.hg38.interval_list}"

SUBMIT="${SUBMIT:-1}"

# --- 覆盖样本清单 / 单样本直投 ---
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

###########################
# 检查
###########################
[[ -z "$ENV_FILE" || -s "$ENV_FILE" ]] || { echo "[FATAL] ENV_FILE not found: $ENV_FILE" >&2; exit 1; }
[[ -s "$REF_FASTA"  ]] || { echo "[FATAL] 缺少参考 fasta: $REF_FASTA" >&2; exit 1; }
[[ -s "$INTERVALS"  ]] || { echo "[FATAL] 缺少 interval_list: $INTERVALS" >&2; exit 1; }
[[ -s "$DBSNP" && -s "${DBSNP}.tbi" ]] || { echo "[FATAL] 缺少 dbSNP 或索引: $DBSNP(.tbi)" >&2; exit 1; }

mkdir -p "${BASE_OUT}/02_haplotypecaller/sh" \
         "${BASE_OUT}/gvcf"

echo "========================================"
echo ">> 胚系 第二步：HaplotypeCaller (GVCF 模式)"
echo "BASE_OUT  : $BASE_OUT"
echo "INTERVALS : $INTERVALS"
echo "QUEUE/CORE/MEM/WALLTIME : $QUEUE / $CORE / ${MEM}G / $WALLTIME"
echo "========================================"

while read -r sample; do
  [[ -z "${sample}" ]] && continue

  shp="${BASE_OUT}/02_haplotypecaller/sh/${sample}.step2.sh"
  cat > "${shp}" <<'EOSH'
#!/usr/bin/env bash
#SBATCH -p __QUEUE__
#SBATCH -J g_step2___SAMPLE__
#SBATCH -c __CORE__
#SBATCH --mem=__MEM__G
#SBATCH -t __WALLTIME__
#SBATCH -o __BASE__/02_haplotypecaller/sh/__SAMPLE__.step2.out
#SBATCH -e __BASE__/02_haplotypecaller/sh/__SAMPLE__.step2.err


set -euo pipefail

[[ -n "__ENV_FILE__" && -s "__ENV_FILE__" ]] && source "__ENV_FILE__"

GATK="__GATK__"
index="__REF_FASTA__"
intervals="__INTERVALS__"
dbsnp="__DBSNP__"

bqsr_bam="__BASE__/BQSR/__SAMPLE__.BQSR.bam"
gvcf="__BASE__/gvcf/__SAMPLE__.g.vcf.gz"

[[ -s "$bqsr_bam" ]] || { echo "缺少 BQSR bam: $bqsr_bam"; exit 1; }

# 已有结果则跳过（避免重复计算）
if [[ -s "$gvcf" && -s "${gvcf}.tbi" ]]; then
  echo "✓ GVCF 已存在，跳过：$gvcf"
  exit 0
fi

echo "[$(date)] __SAMPLE__ 开始 HaplotypeCaller（GVCF 模式）"
"${GATK}" --java-options "-Xms__JAVA_MEM_G__G -Xmx__JAVA_MEM_G__G -XX:ParallelGCThreads=__CORE__" \
  HaplotypeCaller \
  -R "${index}" \
  -L "${intervals}" \
  -I "${bqsr_bam}" \
  -O "${gvcf}" \
  -ERC GVCF \
  --dbsnp "${dbsnp}" \
  --native-pair-hmm-threads __CORE__

echo "[$(date)] 第二步完成：__SAMPLE__ -> ${gvcf}"
EOSH

  sed -i \
    -e "s#__QUEUE__#${QUEUE}#g" \
    -e "s#__SAMPLE__#${sample}#g" \
    -e "s#__CORE__#${CORE}#g" \
    -e "s#__MEM__#${MEM}#g" \
    -e "s#__WALLTIME__#${WALLTIME}#g" \
    -e "s#__JAVA_MEM_G__#${JAVA_MEM_G}#g" \
    -e "s#__BASE__#${BASE_OUT}#g" \
    -e "s#__ENV_FILE__#${ENV_FILE}#g" \
    -e "s#__REF_FASTA__#${REF_FASTA}#g" \
    -e "s#__INTERVALS__#${INTERVALS}#g" \
    -e "s#__DBSNP__#${DBSNP}#g" \
    -e "s#__GATK__#${GATK}#g" \
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
  echo ">> 第二步所有任务已提交完成！"
  echo ">> 等待所有样本 GVCF 生成完毕后，再运行 step3（联合基因分型）"
else
  echo ">> 第二步子脚本已生成完成，但没有提交。"
fi
echo "========================================"
