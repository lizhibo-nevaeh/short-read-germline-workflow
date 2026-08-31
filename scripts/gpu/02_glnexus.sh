#!/usr/bin/env bash
# ============================================================================
# 胚系（Germline）GPU 主流程 第二步
#   工具：GLnexus（DeepVariant GVCF 的联合基因分型工具）
#   作用：把第一步所有样本的 GVCF 合并成一个多样本联合 VCF
#         等价于 CPU 标准流程里的 GenomicsDBImport + GenotypeGVCFs
#
#   注意：
#   - GLnexus 是 CPU 工具（不吃 GPU），不在 Parabricks 镜像内，单独用其 docker 镜像。
#   - 必须使用与 DeepVariant 配套的预设配置：WGS 用 DeepVariantWGS，WES 用 DeepVariantWES。
#   - 本步骤不做 VQSR（DeepVariant 路线不使用 VQSR）；合并后做基础硬过滤
#     （保留 PASS、按 GQ/DP 过滤），具体阈值可按队列调整。
#
#   投递方式（CPU 任务，可在CPU node跑；无调度系统用 nohup）：
#     nohup bash scripts/gpu/02_glnexus.sh /path/to/germline-project/samples.list \
#       > germline_gpu_step2.main.log 2>&1 &
#
#   环境变量覆盖示例：
#     BASE_OUT=/path/to/germline-project SEQ_MODE=wgs CORE=32 MEM=128 \
#     nohup bash scripts/gpu/02_glnexus.sh samples.list > step2.log 2>&1 &
# ============================================================================
set -euo pipefail

###########################
# 配置参数（允许环境变量覆盖）
###########################
SAMPLE_LIST="${1:-${BASE_OUT}/samples.list}"
BASE_OUT="${BASE_OUT:-${PWD}/work}"

# 测序类型决定 GLnexus 预设：wgs -> DeepVariantWGS，wes -> DeepVariantWES
SEQ_MODE="${SEQ_MODE:-wgs}"

# GLnexus docker 镜像
GLNEXUS_IMAGE="${GLNEXUS_IMAGE:-ghcr.io/dnanexus-rnd/glnexus:v1.4.1}"

# bcftools（用于压缩/索引/过滤；从环境取，或用容器内的）
ENV_FILE="${ENV_FILE:-}"
[[ -n "$ENV_FILE" && -s "$ENV_FILE" ]] && source "$ENV_FILE" || true
BCFTOOLS="${BCFTOOLS:-${RESEQ_PREFIX:-/usr/local}/bin/bcftools}"

# 资源
CORE="${CORE:-32}"
MEM="${MEM:-128}"        # GLnexus 内存上限（GB）

# 基础硬过滤阈值（合并后；可按队列调整）
MIN_GQ="${MIN_GQ:-20}"   # 基因型质量
MIN_DP="${MIN_DP:-10}"   # 位点深度

# 输出前缀
COHORT_NAME="${COHORT_NAME:-cohort}"

###########################
# 准备与检查
###########################
sed -i 's/\r$//' "${SAMPLE_LIST}" 2>/dev/null || true
[[ -s "${SAMPLE_LIST}" ]] || { echo "[FATAL] 找不到/为空的样本列表: ${SAMPLE_LIST}" >&2; exit 1; }

case "${SEQ_MODE}" in
  wgs|WGS) GLNEXUS_CONFIG="DeepVariantWGS" ;;
  wes|WES) GLNEXUS_CONFIG="DeepVariantWES" ;;
  *) echo "[FATAL] SEQ_MODE 必须是 wgs 或 wes，当前=${SEQ_MODE}" >&2; exit 1 ;;
esac

mkdir -p "${BASE_OUT}/02_glnexus/log" \
         "${BASE_OUT}"/{joint,final_vcf}

###########################
# 耗时/状态记录（队列级，单个 timing.log）
###########################
TIMING="${BASE_OUT}/02_glnexus/log/${COHORT_NAME}.timing.log"
: > "${TIMING}"
STEP_START_TS=0
fmt_time() { local s=$1; printf '%dh%02dm%02ds' $((s/3600)) $(((s%3600)/60)) $((s%60)); }
step_start() {
  STEP_START_TS=$(date +%s)
  echo "[$(date '+%F %T')] 开始 $1..."
  printf '[%s] START  %s\n' "$(date '+%F %T')" "$1" >> "${TIMING}"
}
step_end() {
  local end_ts elapsed; end_ts=$(date +%s); elapsed=$((end_ts - STEP_START_TS))
  echo "[$(date '+%F %T')] 完成 $1，耗时 $(fmt_time ${elapsed})"
  printf '[%s] END    %-18s elapsed=%-7s (%s)\n' \
    "$(date '+%F %T')" "$1" "${elapsed}s" "$(fmt_time ${elapsed})" >> "${TIMING}"
}

# 收集所有样本 GVCF，生成清单（容器内路径）
gvcf_list_host="${BASE_OUT}/joint/${COHORT_NAME}.gvcf.list"
: > "${gvcf_list_host}"
miss=0; n=0
while read -r sample; do
  [[ -z "${sample}" ]] && continue
  g="${BASE_OUT}/gvcf/${sample}.g.vcf.gz"
  if [[ -s "${g}" ]]; then
    # 写容器内路径（/workdir 映射到 BASE_OUT）
    echo "/workdir/gvcf/${sample}.g.vcf.gz" >> "${gvcf_list_host}"
    n=$((n+1))
  else
    echo "[WARN] 缺少 GVCF: ${g}" >&2
    miss=$((miss+1))
  fi
done < "${SAMPLE_LIST}"

[[ "${n}" -ge 1 ]] || { echo "[FATAL] 没有任何可用 GVCF，请先完成第一步" >&2; exit 1; }
[[ "${miss}" -eq 0 ]] || echo "[WARN] 有 ${miss} 个样本 GVCF 缺失，将仅合并已有的 ${n} 个"

bcf_vcf="${BASE_OUT}/joint/${COHORT_NAME}.glnexus.bcf"
raw_vcf="${BASE_OUT}/joint/${COHORT_NAME}.glnexus.vcf.gz"
pass_vcf="${BASE_OUT}/final_vcf/${COHORT_NAME}.PASS.filtered.vcf.gz"
log_file="${BASE_OUT}/02_glnexus/log/${COHORT_NAME}.log"

echo "=========================================="
echo "胚系 GPU 第二步：GLnexus 联合基因分型"
echo "样本数     : ${n}"
echo "GLnexus配置: ${GLNEXUS_CONFIG}"
echo "镜像       : ${GLNEXUS_IMAGE}"
echo "硬过滤     : GQ>=${MIN_GQ}, DP>=${MIN_DP}"
echo "=========================================="

# GLnexus 要求输出数据库目录不存在，重跑前清理
GLNEXUS_DB="${BASE_OUT}/joint/GLnexus.DB"
case "${GLNEXUS_DB}" in
  "${BASE_OUT}/joint/GLnexus.DB")
    rm -rf -- "${GLNEXUS_DB}"
    ;;
  *)
    echo "[FATAL] Refusing to remove unexpected GLnexus DB path: ${GLNEXUS_DB}" >&2
    exit 1
    ;;
esac

###########################
# 1) GLnexus 合并 -> BCF
###########################
step_start glnexus_merge
docker run --rm \
  -u "$(id -u)":"$(id -g)" \
  --workdir /workdir \
  -v "${BASE_OUT}":/workdir \
  "${GLNEXUS_IMAGE}" /usr/local/bin/glnexus_cli \
    --config "${GLNEXUS_CONFIG}" \
    --dir "/workdir/joint/GLnexus.DB" \
    --threads "${CORE}" \
    --mem-gbytes "${MEM}" \
    --list "/workdir/joint/${COHORT_NAME}.gvcf.list" \
    > "${bcf_vcf}" 2> "${log_file}"

[[ -s "${bcf_vcf}" ]] || { echo "[FATAL] GLnexus 输出为空，详见 ${log_file}" >&2; \
  printf '[%s] FAILED glnexus_merge（输出为空）\n' "$(date '+%F %T')" >> "${TIMING}"; exit 1; }
step_end glnexus_merge

###########################
# 2) BCF -> 压缩 VCF + 索引
###########################
step_start bcf2vcf
"${BCFTOOLS}" view --threads "${CORE}" -Oz -o "${raw_vcf}" "${bcf_vcf}"
"${BCFTOOLS}" index --threads "${CORE}" -t "${raw_vcf}"
step_end bcf2vcf

###########################
# 3) 基础硬过滤（保留 PASS，按 GQ/DP 过滤）
#    DeepVariant 路线不做 VQSR；这里用基因型层面的硬过滤替代
###########################
step_start hard_filter
"${BCFTOOLS}" view --threads "${CORE}" -f PASS "${raw_vcf}" \
  | "${BCFTOOLS}" filter --threads "${CORE}" \
      -e "FMT/GQ<${MIN_GQ} | FMT/DP<${MIN_DP}" \
      --set-GTs . \
      -Oz -o "${pass_vcf}"
"${BCFTOOLS}" index --threads "${CORE}" -t "${pass_vcf}"
step_end hard_filter

printf '[%s] ALL DONE  %s\n' "$(date '+%F %T')" "${COHORT_NAME}" >> "${TIMING}"

echo ""
echo "=========================================="
echo "胚系 GPU 第二步完成！"
echo "  GLnexus 原始联合 VCF : ${raw_vcf}"
echo "  过滤后 PASS VCF      : ${pass_vcf}"
echo "  运行日志             : ${log_file}"
echo "  耗时记录             : ${TIMING}"
echo "  变异数（过滤后）     :"
"${BCFTOOLS}" index -n "${pass_vcf}" 2>/dev/null || zcat "${pass_vcf}" | grep -vc '^#'
echo "=========================================="
