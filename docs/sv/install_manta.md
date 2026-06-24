# Manta Installation and SV Calling

**Reference:** Chen X. et al., Bioinformatics 32(8):1220-1222 (2016)  
**Version used in the manuscript:** Manta 1.6.0

## 1. Install

Manta requires Python 2.7. Keep it in its own Conda environment.

```bash
conda create -n sv_tools -c conda-forge -c bioconda python=2.7 manta samtools
conda activate sv_tools

# Verify
manta --version
```

## 2. Run Manta (Single-Sample Mode)

*Note: NOT used for the manuscript; cannot distinguish radiation-induced from baseline variation.*

```bash
configManta.py \
    --bam        /path/to/HUVEC_DNA_dE_W3.bam \
    --referenceFasta /path/to/GRCh38/Homo_sapiens_assembly38.fasta \
    --runDir     /path/to/output/manta_dE_W3

/path/to/output/manta_dE_W3/runWorkflow.py -m local -j 8
```

## 3. Run Manta (Somatic Mode)

*Note: USED in the manuscript (see Methods §2.4). Treat W0 / d0 as "normal", treated samples as "tumor".*

```bash
configManta.py \
    --normalBam ./bam_files/d0_W1.bam \
    --tumorBam  ./bam_files/dE_W1.bam \
    --referenceFasta /path/to/GRCh38/Homo_sapiens_assembly38.fasta \
    --runDir    ./manta_results/d0_vs_dE_W1

./manta_results/d0_vs_dE_W1/runWorkflow.py -m local -j 8
```

**Output of interest:**
`manta_results/d0_vs_dE_W1/results/variants/somaticSV.vcf.gz`  
This is the file that AnnotSV consumes (see `install_annotsv.md`).

## 4. Batch Run

Batch run for all dose × timepoint combinations.

```bash
REFERENCE=/path/to/GRCh38/Homo_sapiens_assembly38.fasta
OUTBASE=./manta_results

for DOSE in dA dB dC dD dE; do
    for TP in W1 W2 W3; do
        configManta.py \
            --normalBam ./bam_files/d0_${TP}.bam \
            --tumorBam  ./bam_files/${DOSE}_${TP}.bam \
            --referenceFasta ${REFERENCE} \
            --runDir    ${OUTBASE}/d0_vs_${DOSE}_${TP}

        ${OUTBASE}/d0_vs_${DOSE}_${TP}/runWorkflow.py -m local -j 8
    done
done

# Optionally also run d0-vs-d0 controls so the same pipeline can be applied to spontaneous SVs:
for TP in W1 W2 W3; do
    configManta.py \
        --normalBam ./bam_files/d0_W1.bam \
        --tumorBam  ./bam_files/d0_${TP}.bam \
        --referenceFasta ${REFERENCE} \
        --runDir    ${OUTBASE}/d0_vs_d0_${TP}
    ${OUTBASE}/d0_vs_d0_${TP}/runWorkflow.py -m local -j 8
done
```
