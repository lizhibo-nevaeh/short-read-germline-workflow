# Short-Read Germline Workflow

A dual-route workflow for short-read germline small-variant analysis, providing a conventional CPU/GATK route and an accelerated GPU/DeepVariant route.

## Workflow

```text
                         paired-end FASTQ
                               |
               +---------------+---------------+
               |                               |
           CPU route                       GPU route
               |                               |
      fastp + BWA-MEM                 optional fastp
      duplicate marking                      |
            + BQSR                  Parabricks deepvariant_germline
               |                               |
         BQSR BAM                         BAM + GVCF
               |                               |
       HaplotypeCaller                      GLnexus
          GVCF mode                            |
               |                         cohort VCF
      GenomicsDBImport                         |
      + GenotypeGVCFs                 GQ/DP genotype filtering
               |
        chromosome VCFs
               |
          MergeVcfs
               |
       SNP VQSR -> INDEL VQSR
               |
          PASS cohort VCF
```

The two routes are alternatives for generating a cohort-level germline callset. They share the same reference-resource framework but use different variant callers and joint-genotyping strategies.

## Repository structure

```text
short-read-germline-workflow/
├── config/
│   ├── config.example.sh
│   └── samples.example.txt
├── scripts/
│   ├── cpu/
│   │   ├── 01_preprocess_bqsr.sh
│   │   ├── 02_haplotypecaller.sh
│   │   ├── 03_joint_genotyping.sh
│   │   └── 04_vqsr.sh
│   └── gpu/
│       ├── 01_deepvariant.sh
│       └── 02_glnexus.sh
├── .gitignore
└── README.md
```

## Inputs

- paired-end short-read FASTQ files
- one sample ID per line in a sample list
- GRCh38 reference FASTA and indexes
- BQSR known-sites resources
- dbSNP and interval resources
- VQSR training resources for the CPU route

The example scripts expect a sample-oriented layout such as:

```text
data/
├── SAMPLE01/
│   ├── SAMPLE01_f1.fastq.gz
│   └── SAMPLE01_r2.fastq.gz
└── SAMPLE02/
    ├── SAMPLE02_f1.fastq.gz
    └── SAMPLE02_r2.fastq.gz
```

## Configuration

```bash
cp config/config.example.sh config/config.sh
cp config/samples.example.txt /path/to/germline-project/samples.list

# Edit paths and resources, then:
source config/config.sh
```

Software versions and resource choices should be reviewed for the target environment before production use. The container tags in the example configuration are reproducibility examples, not claims of current latest versions.

## CPU route

The CPU route follows a GATK-style cohort workflow.

### 1. Preprocess reads and produce BQSR BAMs

```bash
source config/config.sh
bash scripts/cpu/01_preprocess_bqsr.sh "${SAMPLE_LIST}"
```

Per sample:

```text
fastp -> BWA-MEM -> sort -> duplicate marking -> BaseRecalibrator -> ApplyBQSR
```

The script generates per-sample SLURM jobs by default. Use `SUBMIT=0` to generate job scripts without submitting them.

### 2. HaplotypeCaller in GVCF mode

```bash
source config/config.sh
bash scripts/cpu/02_haplotypecaller.sh "${SAMPLE_LIST}"
```

Produces one GVCF per sample.

### 3. Joint genotyping

```bash
source config/config.sh
bash scripts/cpu/03_joint_genotyping.sh
```

This route imports GVCFs into GenomicsDB by chromosome and runs `GenotypeGVCFs` to generate cohort VCFs.

### 4. Merge and VQSR

```bash
source config/config.sh
bash scripts/cpu/04_vqsr.sh
```

The script merges chromosome-level VCFs, applies SNP VQSR followed by INDEL VQSR, and extracts PASS variants.

## GPU route

The GPU route uses NVIDIA Parabricks DeepVariant followed by GLnexus.

### 1. DeepVariant germline calling

```bash
source config/config.sh
bash scripts/gpu/01_deepvariant.sh "${SAMPLE_LIST}"
```

The supplied script uses `pbrun deepvariant_germline` to produce per-sample BQSR BAMs and GVCFs. Optional fastp preprocessing is controlled by `DO_FASTP`.

### 2. GLnexus joint genotyping

```bash
source config/config.sh
bash scripts/gpu/02_glnexus.sh "${SAMPLE_LIST}"
```

`SEQ_MODE=wgs` selects the `DeepVariantWGS` GLnexus preset and `SEQ_MODE=wes` selects `DeepVariantWES`. The script then applies genotype-level GQ/DP filtering; thresholds can be adjusted for the cohort and validation strategy.

## Main software

- fastp
- BWA-MEM
- samtools
- sambamba
- GATK
- Docker
- NVIDIA Parabricks / DeepVariant
- GLnexus
- bcftools
- SLURM for the CPU route by default

## Notes

- Reference files, sequencing data, software installations, and scheduler profiles are not included.
- Paths, queues, memory, thread counts, interval resources, known-sites resources, and VQSR settings must be adapted to the local environment.
- VQSR requires enough variants/samples for stable model training; small cohorts may require a different filtering strategy.
- The public edition preserves the two route designs and core processing logic while removing environment-specific paths, dataset identifiers, and internal scheduler names.
- Before formal use, validate the chosen software versions, reference resources, thresholds, and a small test dataset in the target environment.
