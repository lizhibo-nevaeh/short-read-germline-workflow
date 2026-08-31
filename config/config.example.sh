#!/usr/bin/env bash
# Example configuration for short-read germline analysis.
# Copy to config/config.sh, edit paths/resources, then source it before running a route.

export BASE_OUT="/path/to/germline-project"
export DATASET_DIR="${BASE_OUT}/data"
export SAMPLE_LIST="${BASE_OUT}/samples.list"

export REF_ROOT="/path/to/reference"
export REF_FASTA="${REF_ROOT}/GRCh38.fasta"
export DBSNP="${REF_ROOT}/dbsnp_146.hg38.vcf.gz"

# BQSR known-sites
export KNOWN1="${REF_ROOT}/dbsnp_146.hg38.vcf.gz"
export KNOWN2="${REF_ROOT}/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
export KNOWN3="${REF_ROOT}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
export KNOWN4="${REF_ROOT}/1000G_omni2.5.hg38.vcf.gz"

# HaplotypeCaller interval list: use a WGS calling-region list or a WES target interval list.
export INTERVALS="${REF_ROOT}/hg38_v0_wgs_calling_regions.hg38.interval_list"

# VQSR resources for the CPU route
export HAPMAP="${REF_ROOT}/hapmap_3.3.hg38.sites.vcf.gz"
export OMNI="${REF_ROOT}/1000G_omni2.5.hg38.vcf.gz"
export SNP_1000G="${REF_ROOT}/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
export MILLS="${REF_ROOT}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
export AXIOM="${REF_ROOT}/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"

# Optional environment bootstrap file. Leave empty when tools are already on PATH.
export ENV_FILE=""
export GATK="gatk"

# Scheduler defaults for the CPU route.
export QUEUE="compute"
export CORE="8"
export MEM="64"
export WALLTIME="24:00:00"
export SUBMIT="1"

# Sequencing mode and read-group library label.
export SEQ_MODE="wgs"   # wgs or wes
export RG_LB="WGS"

# GPU/DeepVariant route
export IMAGE="nvcr.io/nvidia/clara/clara-parabricks:4.5.1-1"
export GLNEXUS_IMAGE="ghcr.io/dnanexus-rnd/glnexus:v1.4.1"

# Joint-calling labels/filter examples.
export COHORT_NAME="example_cohort"
export MIN_GQ="20"
export MIN_DP="10"
