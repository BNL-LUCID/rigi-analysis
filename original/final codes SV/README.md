# SV Analysis Pipeline

Structural-variant half of the LUCID radiation-mutagenesis pipeline.
Takes Manta-called somatic SVs (annotated with AnnotSV) through PASS
filtering, temporal pattern assignment, breakpoint-proximal mutation
enrichment, and dose-stratified functional annotation.

The point-mutation half lives in `../final codes_mutation/` and produces
the per-dose merged CSVs that several scripts here consume.

## Pipeline

```
Manta somaticSV.vcf.gz                    (external — see docs/INSTALL.md)
      ↓
AnnotSV --annotationMode full             (external — see docs/INSTALL.md)
      ↓
01_filter_pass.py             →  PASS-only AnnotSV TSVs
      ↓
02_sv_temporal.py             →  sv_temporal_catalog.csv  (consumed by 5, 6)
      ↓
  ┌───┴───┬──────────────┬──────────────┬───────────────────────┐
  ▼       ▼              ▼              ▼                       ▼
03         04             05             09                      11
SV         repeat         SV-mutation    temporal                dose-stratified
landscape  analysis       correlation    concordance             INV-DBS pairs
(Fig 3)    (Fig S2)       ⭐ 7.13×       (Methods §4.8)          (Fig 5B)
                          (Methods §4.6)        ↓                       ↓
                                ↓               10                      12
                                6               concordance             fetch
                                decay           viz                     annotations
                                (Fig S3)        (Fig S4)                (MyGene+gnomAD)
                                ↓                                       ↓
                                7                                       13
                                INV-size                                categorise
                                ↓                                       genes
                                8                                       ↓
                                INV-DBS figure                          14
                                (Fig 4)                                 dose-response
                                                                        figure (Fig 5)
```

## Setup

```bash
pip install -r requirements.txt
```

External tools (Manta 1.6.0, AnnotSV 3.3) and reference data (gnomAD
v2.1.1) are documented in [`docs/INSTALL.md`](docs/INSTALL.md).

## End-to-end example

Assumes raw AnnotSV outputs are organised into parallel directories
for radiation and control samples (see
[`docs/FILE_LAYOUT.md`](docs/FILE_LAYOUT.md)). Mutation-side merged
CSVs at `../final codes_mutation/merged/`.

```bash
# Step 1: PASS filter (run once per upstream directory)
python 01_filter_pass.py -i ./annotsv         -o ./annotsv_passed
python 01_filter_pass.py -i ./annotsv_control -o ./annotsv_passed_control

# Step 2: Temporal pattern assignment (radiation only)
python 02_sv_temporal.py ./annotsv_passed \
    --output ./out/sv_temporal --tolerance 1000 --plot

# Step 3: SV landscape (Fig 3)
python 03_SV_landscape.py \
    --annotsv-dir ./annotsv_passed \
    --output      ./out/figure3_sv_landscape.png

# Step 4: Repeat analysis (Fig S2) — uses both radiation and control dirs
python 04_repeat_analysis.py \
    --radiation-dir ./annotsv_passed \
    --control-dir   ./annotsv_passed_control \
    --output        ./out/figure_s2_repeat_analysis.png

# Step 5: Breakpoint-proximal mutation enrichment (headline analysis)
python 05_sv_mutation_correlation.py \
    --sv-catalog ./out/sv_temporal/sv_temporal_catalog.csv \
    --mutations  "../final codes_mutation/merged" \
    --output     ./out/sv_correlation \
    --windows    10,25,50,100 --plot

# Step 6: Distance decay (Fig S3)
python 06_sv_type_specific_decay.py \
    --sv-catalog ./out/sv_temporal/sv_temporal_catalog.csv \
    --mutations  "../final codes_mutation/merged/DBS" \
    --output     ./out/sv_type_decay \
    --windows    10,25,50,100

# Step 7: INV-size × DBS coupling (feeds Fig 4 panels C, D)
python 07_inversion_size_analysis.py \
    --annotsv-dir   ./annotsv_passed \
    --dbs-dir       "../final codes_mutation/merged/DBS" \
    --window        10 --mega-threshold 50000000 \
    --output        ./out/inv_size/inv_size.png

# Step 8: INV-DBS coupling figure (Fig 4)
python 08_sv_mut_vizualization.py \
    --correlation-dir ./out/sv_correlation \
    --sv-type-dir     ./out/sv_type_decay \
    --size-analysis   ./out/inv_size \
    --output          ./out/figure4_inv_dbs_unified.png

# Step 9: Temporal concordance (§4.8)
python 09_temporal_concordance.py \
    --annotsv-dir ./annotsv_passed \
    --dbs-data    "../final codes_mutation/merged/DBS" \
    --window      10

# Step 10: Concordance figure (Fig S4)
python 10_temporal_concordance_viz.py \
    --annotsv-dir ./annotsv_passed \
    --dbs-data    "../final codes_mutation/merged/DBS" \
    --output      ./out/figure_temporal_dynamics.png \
    --window      10

# Step 11: Dose-stratified INV-DBS pairs (Fig 5B)
python 11_dose_stratified.py \
    --annotsv-dir    ./annotsv_passed \
    --mutation-dir   "../final codes_mutation/merged" \
    --output-dir     ./out/dose_stratified \
    --window         10 --mega-threshold 50000000

# Step 12: Fetch gene annotations (requires network access)
python 12_fetch_annotations.py \
    --genes-high ./out/dose_stratified/genes_high_dose.csv \
    --genes-low  ./out/dose_stratified/genes_low_dose.csv \
    --gnomad     ./gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz \
    --output     ./out/annotated_genes.csv

# Step 13: Categorise genes (6 functional categories)
python 13_categorise_genes.py \
    --input  ./out/annotated_genes.csv \
    --output ./out/categorized_genes.csv

# Step 14: Dose-response figure (Fig 5)
python 14_dose_based_visualize.py \
    --inv-dbs-low      ./out/dose_stratified/inv_dbs_pairs_low.csv \
    --inv-dbs-high     ./out/dose_stratified/inv_dbs_pairs_high.csv \
    --categorized-genes ./out/categorized_genes.csv \
    --output           ./out/figure5_dose_response
```

## Per-script documentation

Each script's `--help` and module docstring describe its arguments,
inputs, outputs, and edge cases. Highlights:

- **`05_sv_mutation_correlation.py`** is the central analysis (the
  7.13× INV-DBS finding). Runs on Dask-backed dataframes — first
  invocation takes a few minutes.
- **`07_inversion_size_analysis.py`** outputs three CSVs in the
  `--output` parent directory: `size_class_results.csv` (consumed by
  step 8 panel C), `mega_test_results.csv` (Fisher's exact, drives
  panel D significance), and `genes_by_size.csv` (per-gene detail).
- **`08_sv_mut_vizualization.py`** consumes outputs from steps 5, 6,
  and 7 — each panel guards independently, so a missing input directory
  produces "Data not available" rather than failing the whole figure.
- **`12_fetch_annotations.py`** must run with network access (outbound
  HTTPS to MyGene.info).

## Mutation-side prerequisite

Steps 5, 6, 7, 9, 10, and 11 require the merged mutation CSVs from
`../final codes_mutation/`. Run that pipeline first.

## Reference docs

- [`docs/INSTALL.md`](docs/INSTALL.md) — Manta, AnnotSV, gnomAD, MyGene setup
- [`docs/FILE_LAYOUT.md`](docs/FILE_LAYOUT.md) — directory conventions and data flow
- [`docs/FIGURES.md`](docs/FIGURES.md) — manuscript-to-script cross-reference
- [`docs/CONVENTIONS.md`](docs/CONVENTIONS.md) — sample naming, pattern encoding, dose bins

## Citation

If you use this pipeline, please cite the LUCID manuscript and the
upstream tools (Manta, AnnotSV) — full citations in the manuscript
References.
