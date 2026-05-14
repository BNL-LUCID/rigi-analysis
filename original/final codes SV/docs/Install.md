# Install Notes

External tools and reference data that the SV pipeline depends on but
does not bundle.

## Manta 1.6.0

Somatic SV caller. Run upstream of this pipeline; takes BAMs, emits
`somaticSV.vcf.gz`.

Detailed install commands are in `install_manta.txt` (kept as a
standalone file for compatibility with HPC install scripts).

Run mode for this pipeline:

```
configManta.py \
    --normalBam d0_W<n>.bam \
    --tumorBam  d<X>_W<n>.bam \
    --referenceFasta GRCh38.fasta \
    --runDir manta_d<X>_W<n>
manta_d<X>_W<n>/runWorkflow.py
```

Manta must be invoked in **somatic mode** with paired normal/tumor BAMs.
Germline mode produces a different output schema that this pipeline
does not consume.

## AnnotSV 3.3

SV annotator. Runs on Manta's `somaticSV.vcf.gz` and emits the
fully-annotated TSV that all of step 1 onward consumes.

Detailed install commands are in `install_annotsv.txt`.

Run command for this pipeline:

```
AnnotSV \
    -SVinputFile manta_d<X>_W<n>/results/variants/somaticSV.vcf.gz \
    -outputFile  d0_vs_d<X>_W<n>_annotated.tsv \
    -annotationMode full \
    -genomeBuild GRCh38
```

The **`-annotationMode full`** flag is required: without it, AnnotSV
emits one row per gene-overlap, which causes downstream scripts to
multi-count SVs that touch many genes.

## gnomAD constraint file

`12_fetch_annotations.py` parses gnomAD v2.1.1 LoF metrics to assign
pLI scores to genes overlapping concordant INV-DBS pairs.

Download:

```
https://gnomad.broadinstitute.org/downloads
```

File name: `gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz`

Place wherever convenient and pass via `--gnomad`.

## Python dependencies

Python ≥ 3.9. Install via:

```bash
pip install -r requirements.txt
```

Pinned in `requirements.txt`: pandas, numpy, scipy, matplotlib, seaborn,
dask, requests.

## MyGene.info network access

Step 12 (`12_fetch_annotations.py`) makes outbound HTTPS calls to the
MyGene.info REST API for gene name/function lookups. Run outside any
sandboxed or network-restricted environment.
