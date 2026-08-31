#!/usr/bin/env bash
# ========================================
# 胚系流程 第四步 - 合并染色体VCF + VQSR(SNP→InDel串行) + 提取PASS
# ========================================
# 本步骤为“队列级”单次运行（不按样本/染色体拆分）：
#   1) 用 MergeVcfs 把 step3 各染色体的 cohort VCF 合并为整队列 raw VCF
#   2) VQSR 串行：先对 raw VCF 做 SNP 校准并 ApplyVQSR(-mode SNP)，
#      再在其结果上做 InDel 校准并 ApplyVQSR(-mode INDEL) —— GATK Best Practices 推荐方式
#   3) 用 SelectVariants 取 --exclude-filtered 得到 PASS 高置信度变异
#
# 注：默认直接 sbatch 提交一个 step4 作业；也可 SUBMIT=0 只生成子脚本，
#     再用 nohup bash 运行（无 Slurm 服务器适用）。
#
# 投递示例：
# BASE_OUT=/path/to/germline-project \
# REF_ROOT=/path/to/reference \
# ENV_FILE=/path/to/env.sh \
# GATK=gatk \
# CORE=8 MEM=200 SUBMIT=0 \
# bash scripts/cpu/04_vqsr.sh

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

# VQSR 训练资源（示例资源，可按队列与 GATK 建议调整）
HAPMAP="${HAPMAP:-${REF_ROOT}/hapmap_3.3.hg38.sites.vcf.gz}"
OMNI="${OMNI:-${REF_ROOT}/1000G_omni2.5.hg38.vcf.gz}"
SNP_1000G="${SNP_1000G:-${REF_ROOT}/1000G_phase1.snps.high_confidence.hg38.vcf.gz}"
DBSNP="${DBSNP:-${REF_ROOT}/dbsnp_146.hg38.vcf.gz}"
MILLS="${MILLS:-${REF_ROOT}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz}"
AXIOM="${AXIOM:-${REF_ROOT}/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz}"

# VQSR 过滤灵敏度阈值
SNP_TS="${SNP_TS:-99.5}"
INDEL_TS="${INDEL_TS:-99.0}"
INDEL_MAX_GAUSSIANS="${INDEL_MAX_GAUSSIANS:-4}"

# 资源
QUEUE="${QUEUE:-compute}"
CORE="${CORE:-8}"
MEM="${MEM:-200}"
WALLTIME="${WALLTIME:-1000:00:00}"
JAVA_MEM_G="${JAVA_MEM_G:-${MEM}}"

# 染色体列表（需与 step3 一致，用于合并顺序）
CHRS="${CHRS:-chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY}"

SUBMIT="${SUBMIT:-1}"

###########################
# 检查
###########################
[[ -z "$ENV_FILE" || -s "$ENV_FILE" ]] || { echo "[FATAL] ENV_FILE not found: $ENV_FILE" >&2; exit 1; }
[[ -s "$REF_FASTA" ]] || { echo "[FATAL] 缺少参考 fasta: $REF_FASTA" >&2; exit 1; }
for f in "$HAPMAP" "$OMNI" "$SNP_1000G" "$DBSNP" "$MILLS" "$AXIOM"; do
  [[ -s "$f" ]] || { echo "[FATAL] 缺少 VQSR 资源: $f" >&2; exit 1; }
done

mkdir -p "${BASE_OUT}/04_vqsr/sh" \
         "${BASE_OUT}/04_vqsr" \
         "${BASE_OUT}/final_vcf"

# 收集 step3 各染色体 VCF（按 CHRS 顺序），生成 -I 参数列表
INPUT_ARGS=""
miss=0
for chr in ${CHRS}; do
  v="${BASE_OUT}/genotype_vcf/cohort.${chr}.vcf.gz"
  if [[ -s "$v" ]]; then
    INPUT_ARGS="${INPUT_ARGS} -I ${v}"
  else
    echo "[WARN] 缺少染色体 VCF：$v"
    miss=$((miss+1))
  fi
done
[[ -n "$INPUT_ARGS" ]] || { echo "[FATAL] 没有任何可合并的染色体 VCF，请先完成 step3" >&2; exit 1; }
[[ "$miss" -eq 0 ]] || echo "[WARN] 有 $miss 条染色体缺失，将仅合并已有的部分"

echo "========================================"
echo ">> 胚系 第四步：MergeVcfs + VQSR + 提取PASS"
echo "BASE_OUT : $BASE_OUT"
echo "SNP_TS / INDEL_TS : $SNP_TS / $INDEL_TS"
echo "QUEUE/CORE/MEM/WALLTIME : $QUEUE / $CORE / ${MEM}G / $WALLTIME"
echo "========================================"

shp="${BASE_OUT}/04_vqsr/sh/cohort.step4.sh"
cat > "${shp}" <<'EOSH'
#!/usr/bin/env bash
#SBATCH -p __QUEUE__
#SBATCH -J g_step4_vqsr
#SBATCH -c __CORE__
#SBATCH --mem=__MEM__G
#SBATCH -t __WALLTIME__
#SBATCH -o __BASE__/04_vqsr/sh/cohort.step4.out
#SBATCH -e __BASE__/04_vqsr/sh/cohort.step4.err


set -euo pipefail

[[ -n "__ENV_FILE__" && -s "__ENV_FILE__" ]] && source "__ENV_FILE__"

GATK="__GATK__"
index="__REF_FASTA__"

hapmap="__HAPMAP__"
omni="__OMNI__"
snp1000g="__SNP_1000G__"
dbsnp="__DBSNP__"
mills="__MILLS__"
axiom="__AXIOM__"

raw_vcf="__BASE__/04_vqsr/cohort.raw.vcf.gz"
snp_recal="__BASE__/04_vqsr/cohort.snp.recal"
snp_tranches="__BASE__/04_vqsr/cohort.snp.tranches"
snp_rscript="__BASE__/04_vqsr/cohort.snp.plots.R"
snp_filtered="__BASE__/04_vqsr/cohort.snp_filtered.vcf.gz"
indel_recal="__BASE__/04_vqsr/cohort.indel.recal"
indel_tranches="__BASE__/04_vqsr/cohort.indel.tranches"
indel_rscript="__BASE__/04_vqsr/cohort.indel.plots.R"
vqsr_vcf="__BASE__/04_vqsr/cohort.vqsr.vcf.gz"
pass_vcf="__BASE__/final_vcf/cohort.PASS.vcf.gz"

##############################
# 1) 合并各染色体 VCF -> raw VCF
##############################
echo "[$(date)] 开始 MergeVcfs 合并染色体 VCF"
"${GATK}" --java-options "-Xms__JAVA_MEM_G__G -Xmx__JAVA_MEM_G__G -XX:ParallelGCThreads=__CORE__" \
  MergeVcfs \
  __INPUT_ARGS__ \
  -O "${raw_vcf}"

##############################
# 2) SNP 校准模型训练（VariantRecalibrator, -mode SNP）
#    注释：QD MQ DP MQRankSum ReadPosRankSum FS SOR InbreedingCoeff（共8项）
#    InbreedingCoeff 仅适用于样本量≥10的队列
##############################
echo "[$(date)] 开始 SNP VariantRecalibrator"
"${GATK}" --java-options "-Xms__JAVA_MEM_G__G -Xmx__JAVA_MEM_G__G -XX:ParallelGCThreads=__CORE__" \
  VariantRecalibrator \
  -R "${index}" \
  -V "${raw_vcf}" \
  -mode SNP \
  --resource:hapmap,known=false,training=true,truth=true,prior=15.0 "${hapmap}" \
  --resource:omni,known=false,training=true,truth=true,prior=12.0 "${omni}" \
  --resource:1000G,known=false,training=true,truth=false,prior=10.0 "${snp1000g}" \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "${dbsnp}" \
  -an QD -an MQ -an DP -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff \
  -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
  -O "${snp_recal}" \
  --tranches-file "${snp_tranches}" \
  --rscript-file "${snp_rscript}"

echo "[$(date)] 开始 SNP ApplyVQSR (ts=__SNP_TS__)"
"${GATK}" --java-options "-Xms__JAVA_MEM_G__G -Xmx__JAVA_MEM_G__G -XX:ParallelGCThreads=__CORE__" \
  ApplyVQSR \
  -R "${index}" \
  -V "${raw_vcf}" \
  -mode SNP \
  --truth-sensitivity-filter-level __SNP_TS__ \
  --recal-file "${snp_recal}" \
  --tranches-file "${snp_tranches}" \
  -O "${snp_filtered}"

##############################
# 3) InDel 校准模型训练（VariantRecalibrator, -mode INDEL）
#    在 SNP 已校准的 VCF 上串行应用；注释7项 + max-gaussians 4
##############################
echo "[$(date)] 开始 InDel VariantRecalibrator"
"${GATK}" --java-options "-Xms__JAVA_MEM_G__G -Xmx__JAVA_MEM_G__G -XX:ParallelGCThreads=__CORE__" \
  VariantRecalibrator \
  -R "${index}" \
  -V "${snp_filtered}" \
  -mode INDEL \
  --max-gaussians __INDEL_MAX_GAUSSIANS__ \
  --resource:mills,known=false,training=true,truth=true,prior=12.0 "${mills}" \
  --resource:axiomPoly,known=false,training=true,truth=false,prior=10.0 "${axiom}" \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "${dbsnp}" \
  -an QD -an DP -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff \
  -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
  -O "${indel_recal}" \
  --tranches-file "${indel_tranches}" \
  --rscript-file "${indel_rscript}"

echo "[$(date)] 开始 InDel ApplyVQSR (ts=__INDEL_TS__)"
"${GATK}" --java-options "-Xms__JAVA_MEM_G__G -Xmx__JAVA_MEM_G__G -XX:ParallelGCThreads=__CORE__" \
  ApplyVQSR \
  -R "${index}" \
  -V "${snp_filtered}" \
  -mode INDEL \
  --truth-sensitivity-filter-level __INDEL_TS__ \
  --recal-file "${indel_recal}" \
  --tranches-file "${indel_tranches}" \
  -O "${vqsr_vcf}"

##############################
# 4) 提取 PASS 高置信度变异
##############################
echo "[$(date)] 提取 PASS 变异"
"${GATK}" --java-options "-Xms__JAVA_MEM_G__G -Xmx__JAVA_MEM_G__G -XX:ParallelGCThreads=__CORE__" \
  SelectVariants \
  -R "${index}" \
  -V "${vqsr_vcf}" \
  --exclude-filtered \
  -O "${pass_vcf}"

echo "[$(date)] 第四步完成：最终 PASS 变异 -> ${pass_vcf}"
echo "  中间 VQSR 结果（含被过滤记录）：${vqsr_vcf}"
EOSH

sed -i \
  -e "s#__QUEUE__#${QUEUE}#g" \
  -e "s#__CORE__#${CORE}#g" \
  -e "s#__MEM__#${MEM}#g" \
  -e "s#__WALLTIME__#${WALLTIME}#g" \
  -e "s#__JAVA_MEM_G__#${JAVA_MEM_G}#g" \
  -e "s#__BASE__#${BASE_OUT}#g" \
  -e "s#__ENV_FILE__#${ENV_FILE}#g" \
  -e "s#__REF_FASTA__#${REF_FASTA}#g" \
  -e "s#__HAPMAP__#${HAPMAP}#g" \
  -e "s#__OMNI__#${OMNI}#g" \
  -e "s#__SNP_1000G__#${SNP_1000G}#g" \
  -e "s#__DBSNP__#${DBSNP}#g" \
  -e "s#__MILLS__#${MILLS}#g" \
  -e "s#__AXIOM__#${AXIOM}#g" \
  -e "s#__SNP_TS__#${SNP_TS}#g" \
  -e "s#__INDEL_TS__#${INDEL_TS}#g" \
  -e "s#__INDEL_MAX_GAUSSIANS__#${INDEL_MAX_GAUSSIANS}#g" \
  -e "s#__GATK__#${GATK}#g" \
  "${shp}"

# INPUT_ARGS 含斜杠/空格，单独用 # 分隔符替换
sed -i "s#__INPUT_ARGS__#${INPUT_ARGS}#g" "${shp}"

chmod +x "${shp}"
echo "→ 生成 step4 子脚本：$shp"

if [[ "${SUBMIT}" == "1" ]]; then
  if command -v sbatch >/dev/null 2>&1; then
    jid=$(sbatch "${shp}" | awk '{print $4}')
    echo "提交：${shp} -> Job ID: ${jid}"
  else
    echo "[FATAL] SUBMIT=1 但当前环境找不到 sbatch，请用 SUBMIT=0。" >&2
    exit 2
  fi
else
  echo "SUBMIT=0：只生成不提交。可运行："
  echo "  nohup bash ${shp} > ${BASE_OUT}/04_vqsr/sh/cohort.step4.log 2>&1 &"
fi
echo "========================================"
