#!/usr/bin/env bash
# ========================================
# 胚系流程 第三步 - GenomicsDBImport + GenotypeGVCFs（联合基因分型）
# ========================================
# 本步骤为“按染色体”投递（队列级，与 step1/step2 的按样本不同）：
#   1) 自动扫描 gvcf 目录，生成 sample-name-map（样本名 \t GVCF路径）
#   2) 每条染色体生成一个子脚本：GenomicsDBImport 导入该染色体所有样本的 GVCF，
#      再用 GenotypeGVCFs 进行联合基因分型，输出每条染色体的多样本 VCF
#   3) 各染色体跑完后，由 step4 合并并做 VQSR
#
# 投递示例：
# BASE_OUT=/path/to/germline-project \
# REF_ROOT=/path/to/reference \
# ENV_FILE=/path/to/env.sh \
# GATK=gatk \
# CHRS="chr1 chr2 ... chr22 chrX chrY" \
# CORE=8 MEM=200 SUBMIT=0 \
# bash scripts/cpu/03_joint_genotyping.sh

set -euo pipefail
[ -n "${BASH_VERSION:-}" ] || { echo "请用: bash $0"; exit 2; }

###########################
# 配置参数（允许 env 覆盖）
###########################
BASE_OUT="${BASE_OUT:-${PWD}/work}"
ENV_FILE="${ENV_FILE:-}"
REF_ROOT="${REF_ROOT:-${BASE_OUT}/reference}"
GATK="${GATK:-gatk}"

REF_FASTA="${REF_FASTA:-${REF_ROOT}/GRCh38.fasta}"
DBSNP="${DBSNP:-${REF_ROOT}/dbsnp_146.hg38.vcf.gz}"

# 资源（联合基因分型内存需求较高，GenomicsDBImport 建议大内存）
QUEUE="${QUEUE:-compute}"
CORE="${CORE:-8}"
MEM="${MEM:-200}"
WALLTIME="${WALLTIME:-1000:00:00}"
JAVA_MEM_G="${JAVA_MEM_G:-${MEM}}"

# GenomicsDBImport 每批导入样本数 & reader 线程
BATCH_SIZE="${BATCH_SIZE:-50}"
READER_THREADS="${READER_THREADS:-3}"

# GenotypeGVCFs 参数（可按队列调整）
STAND_CALL_CONF="${STAND_CALL_CONF:-30}"
MAX_ALT_ALLELES="${MAX_ALT_ALLELES:-20}"

# 染色体列表（hg38 主染色体；可用 env 覆盖，例如只跑常染色体）
CHRS="${CHRS:-chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY}"

SUBMIT="${SUBMIT:-1}"

# sample-name-map：默认自动生成；如已有可用 SAMPLE_MAP 覆盖
SAMPLE_MAP="${SAMPLE_MAP:-${BASE_OUT}/03_genotype/sample_name_map.txt}"

###########################
# 检查
###########################
[[ -z "$ENV_FILE" || -s "$ENV_FILE" ]] || { echo "[FATAL] ENV_FILE not found: $ENV_FILE" >&2; exit 1; }
[[ -s "$REF_FASTA" ]] || { echo "[FATAL] 缺少参考 fasta: $REF_FASTA" >&2; exit 1; }
[[ -s "$DBSNP" && -s "${DBSNP}.tbi" ]] || { echo "[FATAL] 缺少 dbSNP 或索引: $DBSNP(.tbi)" >&2; exit 1; }

mkdir -p "${BASE_OUT}/03_genotype/sh" \
         "${BASE_OUT}/03_genotype/db" \
         "${BASE_OUT}/genotype_vcf"

###########################
# 生成 sample-name-map（样本名 \t GVCF路径）
###########################
if [[ -s "$SAMPLE_MAP" ]]; then
  echo ">> 使用已存在的 sample-name-map: $SAMPLE_MAP"
else
  echo ">> 扫描 GVCF 目录生成 sample-name-map ..."
  : > "$SAMPLE_MAP"
  shopt -s nullglob
  for g in "${BASE_OUT}/gvcf/"*.g.vcf.gz; do
    s=$(basename "$g" .g.vcf.gz)
    printf '%s\t%s\n' "$s" "$g" >> "$SAMPLE_MAP"
  done
  shopt -u nullglob
fi

n_samples=$(wc -l < "$SAMPLE_MAP" | awk '{print $1}')
[[ "$n_samples" -ge 1 ]] || { echo "[FATAL] sample-name-map 为空，请先完成 step2: $SAMPLE_MAP" >&2; exit 1; }

echo "========================================"
echo ">> 胚系 第三步：GenomicsDBImport + GenotypeGVCFs（联合基因分型）"
echo "样本数      : $n_samples"
echo "染色体      : $CHRS"
echo "BATCH_SIZE  : $BATCH_SIZE"
echo "QUEUE/CORE/MEM/WALLTIME : $QUEUE / $CORE / ${MEM}G / $WALLTIME"
echo "========================================"

gen=0
for chr in ${CHRS}; do
  shp="${BASE_OUT}/03_genotype/sh/${chr}.step3.sh"
  cat > "${shp}" <<'EOSH'
#!/usr/bin/env bash
#SBATCH -p __QUEUE__
#SBATCH -J g_step3___CHR__
#SBATCH -c __CORE__
#SBATCH --mem=__MEM__G
#SBATCH -t __WALLTIME__
#SBATCH -o __BASE__/03_genotype/sh/__CHR__.step3.out
#SBATCH -e __BASE__/03_genotype/sh/__CHR__.step3.err


set -euo pipefail

[[ -n "__ENV_FILE__" && -s "__ENV_FILE__" ]] && source "__ENV_FILE__"

GATK="__GATK__"
index="__REF_FASTA__"
dbsnp="__DBSNP__"
sample_map="__SAMPLE_MAP__"
chr="__CHR__"

db_path="__BASE__/03_genotype/db/${chr}.gdb"
out_vcf="__BASE__/genotype_vcf/cohort.${chr}.vcf.gz"

# 已有结果则跳过
if [[ -s "$out_vcf" && -s "${out_vcf}.tbi" ]]; then
  echo "✓ 该染色体联合分型结果已存在，跳过：$out_vcf"
  exit 0
fi

# GenomicsDBImport 要求工作目录不存在，重跑前先清理
if [[ -d "$db_path" ]]; then
  case "$db_path" in
    "__BASE__/03_genotype/db/"*.gdb)
      echo "[INFO] 清理已存在的 GenomicsDB 目录：$db_path"
      rm -rf -- "$db_path"
      ;;
    *)
      echo "[FATAL] Refusing to remove unexpected GenomicsDB path: $db_path" >&2
      exit 1
      ;;
  esac
fi

echo "[$(date)] ${chr} 开始 GenomicsDBImport"
"${GATK}" --java-options "-Xms__JAVA_MEM_G__G -Xmx__JAVA_MEM_G__G -XX:ParallelGCThreads=__CORE__" \
  GenomicsDBImport \
  --genomicsdb-workspace-path "${db_path}" \
  -L "${chr}" \
  --sample-name-map "${sample_map}" \
  --batch-size __BATCH_SIZE__ \
  --reader-threads __READER_THREADS__

echo "[$(date)] ${chr} 开始 GenotypeGVCFs"
"${GATK}" --java-options "-Xms__JAVA_MEM_G__G -Xmx__JAVA_MEM_G__G -XX:ParallelGCThreads=__CORE__" \
  GenotypeGVCFs \
  -R "${index}" \
  -V "gendb://${db_path}" \
  -O "${out_vcf}" \
  --dbsnp "${dbsnp}" \
  -stand-call-conf __STAND_CALL_CONF__ \
  --max-alternate-alleles __MAX_ALT_ALLELES__

echo "[$(date)] 第三步完成：${chr} -> ${out_vcf}"
EOSH

  sed -i \
    -e "s#__QUEUE__#${QUEUE}#g" \
    -e "s#__CHR__#${chr}#g" \
    -e "s#__CORE__#${CORE}#g" \
    -e "s#__MEM__#${MEM}#g" \
    -e "s#__WALLTIME__#${WALLTIME}#g" \
    -e "s#__JAVA_MEM_G__#${JAVA_MEM_G}#g" \
    -e "s#__BASE__#${BASE_OUT}#g" \
    -e "s#__ENV_FILE__#${ENV_FILE}#g" \
    -e "s#__REF_FASTA__#${REF_FASTA}#g" \
    -e "s#__DBSNP__#${DBSNP}#g" \
    -e "s#__SAMPLE_MAP__#${SAMPLE_MAP}#g" \
    -e "s#__BATCH_SIZE__#${BATCH_SIZE}#g" \
    -e "s#__READER_THREADS__#${READER_THREADS}#g" \
    -e "s#__STAND_CALL_CONF__#${STAND_CALL_CONF}#g" \
    -e "s#__MAX_ALT_ALLELES__#${MAX_ALT_ALLELES}#g" \
    -e "s#__GATK__#${GATK}#g" \
    "${shp}"

  chmod +x "${shp}"
  echo "→ 生成 step3 子脚本：$shp"
  gen=$((gen+1))
done

echo ""
echo "== 生成完成：$gen 条染色体 =="

if [[ "${SUBMIT}" == "1" ]]; then
  if command -v sbatch >/dev/null 2>&1; then
    echo "== 开始批量提交 step3 =="
    sub=0
    for s in "${BASE_OUT}/03_genotype/sh/"*.step3.sh; do
      [[ -s "$s" ]] || continue
      jid=$(sbatch "$s" | awk '{print $4}') || { echo "[WARN] sbatch 失败 -> $s"; continue; }
      echo "提交：$s -> $jid"
      sub=$((sub+1))
    done
    echo "step3 提交完成：$sub / $gen"
    echo ">> 等待所有染色体完成后，再运行 step4（VQSR）"
  else
    echo "[FATAL] SUBMIT=1 但当前环境找不到 sbatch，请用 SUBMIT=0。" >&2
    exit 2
  fi
else
  echo "SUBMIT=0：只生成不提交（可手动 sbatch ${BASE_OUT}/03_genotype/sh/*.step3.sh）"
fi
echo "========================================"
