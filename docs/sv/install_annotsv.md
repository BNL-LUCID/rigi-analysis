# AnnotSV Installation and SV Annotation

**Reference:** Geoffroy V. et al., Bioinformatics 34(20):3572-3574 (2018)
**Version used in the manuscript:** AnnotSV 3.3 (genome build GRCh38)

## 1. Install

Keep AnnotSV in its own Conda environment — its TCL/openblas dependency tends to clash with other genomics environments.

```bash
conda create -n annotsv -c bioconda annotsv
conda activate annotsv
conda install -c conda-forge openblas

# Download the AnnotSV annotation database (GRCh38) following the instructions at
#   https://lbgi.fr/AnnotSV/
# Set the path to the unpacked AnnotSV/ directory below.

ANNOTSV_DB=/path/to/AnnotSV_annotations/share/AnnotSV
```

## 2. Annotate a Single Manta Somatic VCF

`-annotationMode full` → one fully-annotated row per SV (REQUIRED by all downstream scripts in this pipeline; without it, INVs are reported with one row per affected gene, inflating counts 20-100×).

```bash
AnnotSV \
    -SVinputFile     ./manta_results/d0_vs_dE_W1/results/variants/somaticSV.vcf.gz \
    -annotationsDir  ${ANNOTSV_DB} \
    -outputFile      d0_vs_dE_W1_annotated \
    -genomeBuild     GRCh38 \
    -SVinputInfo     1 \
    -annotationMode  full
```

## 3. Batch Annotate Every Manta Somatic VCF

```bash
INPUT_BASE=./manta_results
OUTPUT_DIR=./SV_files/annotated
mkdir -p ${OUTPUT_DIR}

# Radiation samples (d0 vs dA-dE × W1-W3)
for DOSE in dA dB dC dD dE; do
    for TP in W1 W2 W3; do
        AnnotSV \
            -SVinputFile     ${INPUT_BASE}/d0_vs_${DOSE}_${TP}/results/variants/somaticSV.vcf.gz \
            -annotationsDir  ${ANNOTSV_DB} \
            -outputFile      ${OUTPUT_DIR}/d0_vs_${DOSE}_${TP}_annotated \
            -genomeBuild     GRCh38 \
            -SVinputInfo     1 \
            -annotationMode  full
    done
done

# Control samples (d0 vs d0 × W1-W3) — only needed if you want spontaneous SVs
for TP in W1 W2 W3; do
    AnnotSV \
        -SVinputFile     ${INPUT_BASE}/d0_vs_d0_${TP}/results/variants/somaticSV.vcf.gz \
        -annotationsDir  ${ANNOTSV_DB} \
        -outputFile      ${OUTPUT_DIR}/d0_vs_d0_${TP}_annotated \
        -genomeBuild     GRCh38 \
        -SVinputInfo     1 \
        -annotationMode  full
    done
```

## 4. Output

Output is one TSV per Manta run, named e.g. `d0_vs_dE_W1_annotated.tsv`.

These TSVs are the input directory expected by `filter_pass` (see [README.md](README.md) for the downstream pipeline).
