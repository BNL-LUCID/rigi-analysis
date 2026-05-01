# SV Analysis Pipeline

Structural-variant analysis pipeline for the manuscript
**"A Spatio-Temporal Analysis Framework for Characterizing Radiation-Induced Genomic Instability"** (Chopra et al.).

This pipeline takes Manta-called somatic SVs (annotated with AnnotSV) and
produces the temporal pattern catalogues, breakpoint-proximal mutation
enrichment statistics, dose-stratified analyses, and figures presented in the
manuscript (Methods §§2.4–2.9, Figures 3–5, Supplementary Figures S2–S4).

The point-mutation side of the pipeline lives in `../final codes_mutation/` and
produces the per-mutation `*_dose_*_merged.csv` files that the breakpoint-proximal
analyses here consume.

---

## Pipeline overview

```
            Manta somaticSV.vcf.gz                  (external — see install_manta.txt)
                       │
                       ▼
           AnnotSV --annotationMode full            (external — see install_annotsv.txt)
                       │
                       ▼
    01_filter_pass.py     →  PASS-only AnnotSV TSVs
                       │
                       ▼
    02_sv_temporal.py     →  sv_temporal_catalog.csv      (Methods §2.5, Fig S4)
                       │
        ┌──────────────┼──────────────┬──────────────────────────────┐
        ▼              ▼              ▼                              ▼
  03_SV_landscape    04_repeat   05_sv_mutation_correlation     09_temporal_concordance
   (Fig 3)          _analysis    (Methods §2.6,   ⭐ 7.13×)       (Methods §2.8)
                    (Fig S2)            │                              │
                                        ▼                              ▼
                              06_sv_type_specific_decay   10_temporal_concordance_viz
                                  (Fig S3)                 (Fig S4)
                                        │
                                        ▼
                              07_inversion_size_analysis  (INV-size × DBS coupling,
                                                           feeds Fig 4 panels C, D)
                                        │
                                        ▼
                              08_sv_mut_vizualization
                                  (Fig 4, consumes 05, 06 and 07)

                                        │
                                        ▼
                            11_dose_stratified            (Methods §2.9, Fig 5B)
                                        │
                                        ▼
                            12_fetch_annotations          (MyGene + gnomAD pLI)
                                        │
                                        ▼
                            13_categorise_genes           (6 functional categories)
                                        │
                                        ▼
                            14_dose_based_visualize       (Fig 5A, 5B, 5C)
```

`fast_dbs_search.py` is a shared helper module (parallel binary-search DBS
counter) that downstream scripts can import.

---

## Files

| File | Purpose | Manuscript section / figure |
|------|---------|-----------------------------|
| `install_manta.txt`              | Manta install + somatic-mode SV calling commands | §2.4 |
| `install_annotsv.txt`            | AnnotSV install + `-annotationMode full` commands | §2.4 |
| `fast_dbs_search.py`          | Shared parallel DBS-near-breakpoint counter (importable helper) | — |
| `01_filter_pass.py`              | Retain MANTA `FILTER == PASS` rows in every AnnotSV TSV | §2.4 |
| `02_sv_temporal.py`              | 4-position temporal pattern (W0/W1/W2/W3, ±1000 bp) → `sv_temporal_catalog.csv` | §2.5, Fig S4 |
| `03_SV_landscape.py`             | SV counts by type, dose, timepoint; dose × timepoint heatmap | Fig 3 |
| `04_repeat_analysis.py`          | Repeat-element involvement at SV breakpoints | §2.7, Fig S2 |
| `05_sv_mutation_correlation.py`  | **Central analysis.** Mutation enrichment in 10/25/50/100 bp windows around all SV breakpoints (4 SV types × 4 mutation types) | §2.6, Fig 4 |
| `06_sv_type_specific_decay.py`   | Distance-decay curves per SV type | Fig S3 |
| `07_inversion_size_analysis.py`  | INV-size-class × DBS coupling (Tiny → Mega ≥50 Mb), Fisher's test for Mega vs smaller; feeds Fig 4 panels C and D | §2.6, Fig 4 |
| `08_sv_mut_vizualization.py`     | INV–DBS coupling figure (4 panels of evidence; consumes outputs of 05, 06, 07) | Fig 4 |
| `09_temporal_concordance.py`     | INV-timepoint vs DBS-pattern concordance (33 % null) | §2.8 |
| `10_temporal_concordance_viz.py` | Concordance figure | Fig S4 |
| `11_dose_stratified.py`          | Low-bin (A,B,C) vs High-bin (D,E) INV-DBS pairs requiring spatial + dose + temporal concordance | §2.9, Fig 5B |
| `12_fetch_annotations.py`        | MyGene.info REST queries + gnomAD pLI parsing for affected genes | §2.9 |
| `13_categorise_genes.py`         | 6-category functional classification with priority hierarchy | §2.9, Fig 5C |
| `14_dose_based_visualize.py`     | Dose-response figure (genome-wide distribution, INV-DBS counts, functional categories) | Fig 5 |
| `fig1.md`                        | Mermaid source for Figure 1 (workflow schematic) | Fig 1 |

---

## Prerequisites

### External tools (not bundled)

| Tool      | Version | Purpose                                  | Install instructions    |
|-----------|---------|------------------------------------------|-------------------------|
| Manta     | 1.6.0   | Somatic SV calling from BAMs             | `install_manta.txt`     |
| AnnotSV   | 3.3     | SV annotation, repeat overlap            | `install_annotsv.txt`   |
| samtools  | ≥1.10   | BAM handling for Manta                   | (bundled with Manta)    |

Manta must be run in **somatic mode** (`--normalBam d0_Wn.bam --tumorBam dX_Wn.bam`)
and AnnotSV must be invoked with `-annotationMode full` so each SV produces a
single fully-annotated row. See the manuscript Methods §2.4 for parameter
details.

### External data (downloaded locally)

`12_fetch_annotations.py` requires the gnomAD v2.1.1 constraint file:

```
gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz
```

Download from <https://gnomad.broadinstitute.org/downloads>.

### Python

Python ≥ 3.9. Install Python deps with:

```bash
pip install -r requirements.txt
```

### Mutation-side inputs

Several scripts (05, 06, 07, 09, 10, 11, figures 4 & 5) take DBS / SNV / MNS / ID
mutation tables produced by the sister pipeline at `../final codes_mutation/`.
Run that pipeline first to generate the `*_dose_*_merged.csv` files.

---

## Running the pipeline

The example commands assume the AnnotSV outputs are organised into two
parallel directories — one for radiation samples, one for unirradiated
controls:

```
annotsv/                    raw radiation TSVs (d0_vs_d{A,B,C,D,E}_W{1,2,3}_annotated.tsv)
annotsv_control/            raw control TSVs   (d0_vs_d0_W{1,2,3}_annotated.tsv)
```

After `01_filter_pass.py`, the corresponding PASS-only outputs are:

```
annotsv_passed/             PASS-filtered radiation TSVs
annotsv_passed_control/     PASS-filtered control TSVs
```

Mutation merged CSVs (used from step 5 onwards) are expected at
`../final codes_mutation/merged/`.

### Step 1 — PASS filter (run twice, once per upstream directory)

```bash
# Radiation samples
python 01_filter_pass.py \
    --input-dir   ./annotsv \
    --output-dir  ./annotsv_passed

# Control samples
python 01_filter_pass.py \
    --input-dir   ./annotsv_control \
    --output-dir  ./annotsv_passed_control
```

All downstream scripts run on the PASS-filtered directories.

### Step 2 — Temporal pattern assignment (§2.5)

Only uses the radiation directory (it iterates through doses A–E by
construction; controls are handled separately in step 4):

```bash
python 02_sv_temporal.py ./annotsv_passed \
    --output    ./out/sv_temporal \
    --tolerance 1000 \
    --plot
```

Produces `sv_temporal_catalog.csv` (input for many downstream steps), the
per-dose summary, and the temporal-pattern plot. `02_sv_temporal.py` filters
to `Annotation_mode == 'full'` internally. Default tolerance is 1000 bp,
matching manuscript §2.5; `--tolerance 1000` is the manuscript value.

### Step 3 — SV landscape (Fig 3)

Radiation directory only:

```bash
python 03_SV_landscape.py \
    --annotsv-dir ./annotsv_passed \
    --output      ./out/figure3_sv_landscape.png
```

### Step 4 — Repeat element analysis (§2.7, Fig S2)

This is the only step in the upstream-only block that needs **both**
directories (Panel D compares radiation vs control repeat involvement):

```bash
python 04_repeat_analysis.py \
    --radiation-dir ./annotsv_passed \
    --control-dir   ./annotsv_passed_control \
    --output        ./out/figure_s2_repeat_analysis.png
```

> **Known limitations of the bundled `04_repeat_analysis.py`** (worth a once-over before
> trusting the numbers against the manuscript):
>
> 1. **Panels C and D count `Has_Left_Repeat` only**, not "either breakpoint".
>    The manuscript reports per-SV-type and overall repeat % using "any
>    breakpoint in a repeat" semantics (e.g. radiation 73.7%, control 77.9%
>    overall; INV 84.8% per type). Using left only systematically undercounts
>    by a few percent.
> 2. **Panel A pie-chart colors are mismatched** — the colors `[Other, SINE, None]`
>    are repeat-class colors used as labels for the involvement classes
>    `Both in Repeat / One in Repeat / No Repeat`. Visually OK, semantically wrong.
> 3. **No statistical test.** The manuscript reports Fisher's exact test
>    (OR=1.26, p<0.001) for the radiation-vs-control comparison; the script
>    doesn't compute it.
>
> All three are small fixes — flag if you want them done before figure
> finalisation.

### Step 5 — Breakpoint-proximal mutation enrichment (§2.6, central analysis)

```bash
python 05_sv_mutation_correlation.py \
    --sv-catalog ./out/sv_temporal/sv_temporal_catalog.csv \
    --mutations  "../final codes_mutation/merged" \
    --output     ./out/sv_correlation \
    --windows    10,25,50,100 \
    --plot
```

This is the analysis that produces the headline 7.13× INV-DBS enrichment at
10 bp. It runs on Dask-backed dataframes — be patient on the first invocation.

### Step 6 — SV-type-specific distance decay (Fig S3)

```bash
python 06_sv_type_specific_decay.py \
    --sv-catalog ./out/sv_temporal/sv_temporal_catalog.csv \
    --mutations  "../final codes_mutation/merged/DBS" \
    --output     ./out/sv_type_decay \
    --windows    10,25,50,100
```

### Step 7 — Inversion size-class × DBS coupling (feeds Fig 4 panels C, D)

Quantifies the fraction of genes within each INV size class (Tiny < 1 kb,
Small 1-10 kb, Medium 10-100 kb, Large 0.1-1 Mb, Very Large 1-10 Mb,
Huge 10-50 Mb, Mega ≥50 Mb) that contain a DBS within `--window` bp of an
INV breakpoint, plus a Fisher's exact test of Mega vs <Mega.

```bash
python 07_inversion_size_analysis.py \
    --annotsv-dir   ./annotsv_passed \
    --dbs-dir       "../final codes_mutation/merged/DBS" \
    --window        10 \
    --mega-threshold 50000000 \
    --output        ./out/inv_size/inv_size.png
```

Outputs (in the `--output` parent dir, here `./out/inv_size/`):

- `size_class_results.csv` — one row per size class, with `Size_Class`,
  `N_Genes`, `N_With_DBS`, `Percent_With_DBS`. **Step 8 reads this.**
- `mega_test_results.csv` — Fisher's exact for Mega vs <Mega; `P_Value` column
  drives the significance annotation in Step 8 panel D.
- `genes_by_size.csv` — per-gene detail (gene name, max enclosing INV size,
  DBS hit count). Not consumed downstream; kept for inspection.
- `inv_size.png` — standalone size-class bar plot (script's own figure).

> **Pattern filter caveat.** The script counts **all** DBS within the
> window, not just radiation-pattern DBS. If figure 4 panels C/D should be
> radiation-specific, pre-filter `merged/DBS/` to radiation-pattern rows
> before pointing `--dbs-dir` at it, or patch `load_dbs_mutations` to drop
> control-pattern rows. The numbers we ran during validation matched the
> manuscript without filtering — confirm before relying on this default.

### Step 8 — INV-DBS coupling figure (Fig 4)

```bash
python 08_sv_mut_vizualization.py \
    --correlation-dir ./out/sv_correlation \
    --sv-type-dir     ./out/sv_type_decay \
    --size-analysis   ./out/inv_size \
    --output          ./out/figure4_inv_dbs_unified.png
```

Consumes outputs from Steps 5 (panel A), 6 (panel B), and 7 (panels C, D).
Each panel guards independently — a missing input directory leaves that
panel as "Data not available" rather than failing the whole figure.

### Step 9 — Temporal concordance (§2.8)

```bash
python 09_temporal_concordance.py \
    --annotsv-dir ./annotsv_passed \
    --dbs-data    "../final codes_mutation/merged/DBS" \
    --window      10
```

Reports the fraction of co-occurring INV-DBS pairs whose timepoints actually
match (manuscript reports W1 = 66.7 %, W2 = 85.0 %, W3 = 100 % vs the 33 %
null expectation). Reports both "loose" (timepoint anywhere in the DBS
pattern) and "strict" (single-timepoint exact match) definitions.

### Step 10 — Concordance figure (Fig S4)

```bash
python 10_temporal_concordance_viz.py \
    --annotsv-dir ./annotsv_passed \
    --dbs-data    "../final codes_mutation/merged/DBS" \
    --output      ./out/figure_temporal_dynamics.png \
    --window      10
```

### Step 11 — Dose-stratified INV-DBS pairs (§2.9, Fig 5B)

```bash
python 11_dose_stratified.py \
    --annotsv-dir   ./annotsv_passed \
    --mutation-dir  "../final codes_mutation/merged" \
    --output-dir    ./out/dose_stratified \
    --inv-size      0 \
    --window        10 \
    --mega-threshold 50000000
```

Splits doses into Low (A, B, C: 0.001–0.1 mGy/hr) and High (D, E: 1–2 mGy/hr)
bins, requires spatial + dose + temporal concordance, and produces the per-bin
INV-DBS catalogues that feed the functional annotation step.

### Step 12 — Fetch gene annotations (MyGene.info + gnomAD)

```bash
python 12_fetch_annotations.py \
    --genes-high ./out/dose_stratified/genes_high_dose.csv \
    --genes-low  ./out/dose_stratified/genes_low_dose.csv \
    --gnomad     ./gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz \
    --output     ./out/annotated_genes.csv
```

Run this **outside** any sandboxed environment (the script makes outbound
HTTPS calls to the MyGene.info REST API).

### Step 13 — Categorise genes (6 functional categories)

```bash
python 13_categorise_genes.py \
    --input  ./out/annotated_genes.csv \
    --output ./out/categorized_genes.csv
```

Applies the category priority hierarchy from manuscript §2.9: Cell Cycle &
DNA Damage > Signal Transduction > Gene Expression > Cell Structure &
Adhesion > Development & Differentiation > Metabolism & Other.

### Step 14 — Dose-response figure (Fig 5)

```bash
python 14_dose_based_visualize.py \
    --inv-dbs-low      ./out/dose_stratified/inv_dbs_pairs_low.csv \
    --inv-dbs-high     ./out/dose_stratified/inv_dbs_pairs_high.csv \
    --categorized-genes ./out/categorized_genes.csv \
    --output           ./out/figure5_dose_response
```

---

## Method-to-script cross-reference

| Manuscript section / figure | Script(s) |
|---|---|
| §2.4 PASS filter | `01_filter_pass.py` |
| §2.5 Temporal pattern (4-pos C/T/0, ±1000 bp) | `02_sv_temporal.py` |
| §2.6 SV-mutation co-occurrence (Poisson, 10/25/50/100 bp) | `05_sv_mutation_correlation.py`, `06_sv_type_specific_decay.py`, `07_inversion_size_analysis.py` |
| §2.7 Repeat-element analysis | `04_repeat_analysis.py` |
| §2.8 Temporal concordance (33 % null, W1/W2/W3, size strata) | `09_temporal_concordance.py` |
| §2.9 Dose-stratified (Low A-C / High D-E) + functional annotation | `11_dose_stratified.py`, `12_fetch_annotations.py`, `13_categorise_genes.py` |
| Fig 3 SV landscape | `03_SV_landscape.py` |
| Fig 4 INV-DBS coupling | `08_sv_mut_vizualization.py` (consumes 05, 06 & 07 outputs) |
| Fig 5 Dose response | `14_dose_based_visualize.py` (consumes 11 & 13 outputs) |
| Fig S2 Repeat involvement | `04_repeat_analysis.py` |
| Fig S3 Distance decay | `06_sv_type_specific_decay.py` |
| Fig S4 Temporal concordance | `10_temporal_concordance_viz.py` |
| Fig 1 Workflow schematic | `fig1.md` (Mermaid source) |

---

## Conventions

- **Sample naming.** `d0` is unirradiated control. Doses are `dA`–`dE`
  (0.001 / 0.01 / 0.1 / 1 / 2 mGy/hr). Timepoints are `W1`, `W2`, `W3`
  (W0 = pre-exposure baseline, treated as "absent" in the leading position of
  every temporal pattern).
- **Pattern strings.** Four characters, position = `[W0][W1][W2][W3]`.
  Treatment presence = `T`, control presence = `C`, both = `B`, absent = `0`.
  The seven radiation-specific patterns are
  `0T00, 00T0, 000T, 0TT0, 0T0T, 00TT, 0TTT`.
- **AnnotSV `Annotation_mode`.** All scripts that count SVs filter to
  `Annotation_mode == 'full'` so a single SV touching N genes is counted once,
  not N times.
- **Inversion deduplication.** AnnotSV can report the same INV with reversed
  coordinates; the relevant scripts collapse these using sorted `(start, end)`
  keys before counting.
- **Dose bins.** Low = `A, B, C` (0.001–0.1 mGy/hr); High = `D, E` (1–2 mGy/hr).

---

## Citation

Cite the Manta and AnnotSV papers if you use this pipeline:

- Chen X. *et al.* Manta: rapid detection of structural variants and indels
  for germline and cancer sequencing applications. *Bioinformatics* 32,
  1220–1222 (2016).
- Geoffroy V. *et al.* AnnotSV: an integrated tool for structural variations
  annotation. *Bioinformatics* 34, 3572–3574 (2018).

For the analysis framework itself, cite the LUCID manuscript.
