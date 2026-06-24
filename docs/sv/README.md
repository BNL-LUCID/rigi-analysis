# SV Analysis Pipeline (`rigi_analysis.sv`)

Structural-variant half of the LUCID radiation-mutagenesis pipeline.
Takes Manta-called somatic SVs (annotated with AnnotSV) through PASS
filtering, temporal pattern assignment, breakpoint-proximal mutation
enrichment, and dose-stratified functional annotation.

The point-mutation half lives in `src/rigi_analysis/mutation/` and produces
the per-dose merged CSVs that several scripts here consume.

## Pipeline

```
Manta somaticSV.vcf.gz                    (external — see install.md)
      ↓
AnnotSV --annotationMode full             (external — see install.md)
      ↓
filter_pass             →  PASS-only AnnotSV TSVs
      ↓
sv_temporal             →  sv_temporal_catalog.csv  (consumed by 5, 6)
      ↓
  ┌───┴───┬──────────────┬──────────────┬───────────────────────┐
  ▼       ▼              ▼              ▼                       ▼
SV         repeat         SV-mutation    temporal                dose-stratified
landscape  analysis       correlation    concordance             INV-DBS pairs
(Fig 3)    (Fig S2)       ⭐ 7.13×       (Methods §4.8)          (Fig 5B)
                          (Methods §4.6)        ↓                       ↓
                                ↓               concordance             fetch
                                decay           viz                     annotations
                                (Fig S3)        (Fig S4)                (MyGene+gnomAD)
                                ↓                                       ↓
                                INV-size                                categorise
                                ↓                                       genes
                                                                        ↓
                                INV-DBS figure                          dose-response
                                (Fig 4)                                 figure (Fig 5)
```

## Setup

### 1. `rigi-analysis` Environment
The core Python dependencies are managed through a unified Conda environment. Please refer to the **[Main Setup Instructions](../README.md#setup)** for details on creating and activating the `rigi_analysis` environment.

### 2. Third-Party Tools & Data Environment
External tools (Manta 1.6.0, AnnotSV 3.3) and reference data (gnomAD
v2.1.1) are documented in [`install.md`](install.md).

## End-to-end example

Assumes raw AnnotSV outputs are organised into parallel directories
for radiation and control samples. Mutation-side merged
CSVs at `../out_mutation/merged_data/` or similar.

```bash
# Step 1: PASS filter (run once per upstream directory)
rigi-analysis-run filter_pass -i ./annotsv         -o ./annotsv_passed
rigi-analysis-run filter_pass -i ./annotsv_control -o ./annotsv_passed_control

# Step 2: Temporal pattern assignment (radiation only)
rigi-analysis-run sv_temporal ./annotsv_passed \
    --output ./out/sv_temporal --tolerance 1000 --plot

# Step 3: SV landscape (Fig 3)
rigi-analysis-run SV_landscape \
    --annotsv-dir ./annotsv_passed \
    --output      ./out/figure3_sv_landscape.png

# Step 4: Repeat analysis (Fig S2) — uses both radiation and control dirs
rigi-analysis-run repeat_analysis \
    --radiation-dir ./annotsv_passed \
    --control-dir   ./annotsv_passed_control \
    --output        ./out/figure_s2_repeat_analysis.png

# Step 5: Breakpoint-proximal mutation enrichment (headline analysis)
rigi-analysis-run sv_mutation_correlation \
    --sv-catalog ./out/sv_temporal/sv_temporal_catalog.csv \
    --mutations  "../mutation/merged" \
    --output     ./out/sv_correlation \
    --windows    10,25,50,100 --plot

# Step 6: Distance decay (Fig S3)
rigi-analysis-run sv_type_specific_decay \
    --sv-catalog ./out/sv_temporal/sv_temporal_catalog.csv \
    --mutations  "../mutation/merged/DBS" \
    --output     ./out/sv_type_decay \
    --windows    10,25,50,100

# Step 7: INV-size × DBS coupling (feeds Fig 4 panels C, D)
rigi-analysis-run inversion_size_analysis \
    --annotsv-dir   ./annotsv_passed \
    --dbs-dir       "../mutation/merged/DBS" \
    --window        10 --mega-threshold 50000000 \
    --output        ./out/inv_size/inv_size.png

# Step 8: INV-DBS coupling figure (Fig 4)
rigi-analysis-run sv_mut_vizualization \
    --correlation-dir ./out/sv_correlation \
    --sv-type-dir     ./out/sv_type_decay \
    --size-analysis   ./out/inv_size \
    --output          ./out/figure4_inv_dbs_unified.png

# Step 9: Temporal concordance (§4.8)
rigi-analysis-run temporal_concordance \
    --annotsv-dir ./annotsv_passed \
    --dbs-data    "../mutation/merged/DBS" \
    --window      10

# Step 10: Concordance figure (Fig S4)
rigi-analysis-run temporal_concordance_viz \
    --annotsv-dir ./annotsv_passed \
    --dbs-data    "../mutation/merged/DBS" \
    --output      ./out/figure_temporal_dynamics.png \
    --window      10

# Step 11: Dose-stratified INV-DBS pairs (Fig 5B)
rigi-analysis-run dose_stratified \
    --annotsv-dir    ./annotsv_passed \
    --mutation-dir   "../mutation/merged" \
    --output-dir     ./out/dose_stratified \
    --window         10 --mega-threshold 50000000

# Step 12: Fetch gene annotations (requires network access)
rigi-analysis-run fetch_annotations \
    --genes-high ./out/dose_stratified/genes_high_dose.csv \
    --genes-low  ./out/dose_stratified/genes_low_dose.csv \
    --gnomad     ./gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz \
    --output     ./out/annotated_genes.csv

# Step 13: Categorise genes (6 functional categories)
rigi-analysis-run categorise_genes \
    --input  ./out/annotated_genes.csv \
    --output ./out/categorized_genes.csv

# Step 14: Dose-response figure (Fig 5)
rigi-analysis-run dose_based_visualize \
    --inv-dbs-low      ./out/dose_stratified/inv_dbs_pairs_low.csv \
    --inv-dbs-high     ./out/dose_stratified/inv_dbs_pairs_high.csv \
    --categorized-genes ./out/categorized_genes.csv \
    --output           ./out/figure5_dose_response
```

## Per-script documentation

Each script's `--help` and module docstring describe its arguments,
inputs, outputs, and edge cases. Highlights:

- **`sv_mutation_correlation`** is the central analysis (the
  7.13× INV-DBS finding). Runs on Dask-backed dataframes — first
  invocation takes a few minutes.
- **`inversion_size_analysis`** outputs three CSVs in the
  `--output` parent directory: `size_class_results.csv` (consumed by
  step 8 panel C), `mega_test_results.csv` (Fisher's exact, drives
  panel D significance), and `genes_by_size.csv` (per-gene detail).
- **`sv_mut_vizualization`** consumes outputs from steps 5, 6,
  and 7 — each panel guards independently, so a missing input directory
  produces "Data not available" rather than failing the whole figure.
- **`fetch_annotations`** must run with network access (outbound
  HTTPS to MyGene.info).

## Mutation-side prerequisite

Steps 5, 6, 7, 9, 10, and 11 require the merged mutation CSVs from
the mutation pipeline. Run that pipeline first.

## Workflow orchestration

The `SVPipeline` class in `rigi_analysis.workflows.sv_workflow` runs all 14 steps above as an async DAG via RADICAL asyncflow. Steps 1–2 run in parallel; landscape and repeat analysis are independent of the main enrichment chain. To run programmatically:

```python
import asyncio
from rigi_analysis.workflows import SVPipeline

config = {...}  # see examples/ for a full config template
asyncio.run(SVPipeline.run_workflow(config=config))
```

For end-to-end execution of both pipelines in sequence, use `FullWorkflow` in `rigi_analysis.workflows.full_workflow`.

## Reference docs

- [`install.md`](install.md) — Manta, AnnotSV, gnomAD, MyGene setup
- [`conventions.md`](conventions.md) — sample naming, pattern encoding, dose bins
- [`pipeline_diagram.md`](pipeline_diagram.md) — full Mermaid DAG diagram

## Citation

If you use this pipeline, please cite the LUCID manuscript and the
upstream tools (Manta, AnnotSV) — full citations in the manuscript
References.
