# LUCID — Spatio-Temporal Analysis of Radiation-Induced Genomic Instability

Analytical code for the manuscript

> **A Spatio-Temporal Analysis Framework for Characterizing Radiation-Induced Genomic Instability**
> Chopra K., Cucinell C., Weinberg R., Forrester S., Brettin T., Kilic O. O., Yoon B. J.

The pipeline integrates point-mutation profiling and structural-variant
analysis on whole-genome sequencing of HUVEC cells exposed to chronic low-dose
gamma radiation (0.001–2 mGy/hr) over three weeks. It produces the temporal
pattern catalogues, breakpoint-proximal mutation enrichment statistics,
dose-stratified analyses, and figures presented in the paper.

---

## What's in this repository

The pipeline is split into two halves that communicate through one
intermediate artefact:

```
                       Whole-genome sequencing
                       (HUVEC, 5 doses × 3 weekly timepoints)
                                  │
              ┌───────────────────┴───────────────────┐
              ▼                                       ▼
     POINT MUTATIONS                          STRUCTURAL VARIANTS
     Mutect2 (external)                       Manta (external)
              │                                       │
              ▼                                       ▼
     final codes_mutation/                    final codes SV/
   ┌──────────────────────┐                ┌──────────────────────┐
   │ SigProfiler          │                │ AnnotSV --full       │
   │ preprocessing        │                │ PASS filter          │
   │ annotation (refGene) │                │ temporal patterns    │
   │ temporal patterns    │                │ SV landscape (Fig 3) │
   │ Sankey (Fig 2, S1)   │                │ repeat analysis (S2) │
   │ merged CSVs ─────────┼──┐             │                      │
   └──────────────────────┘  │             │                      │
                             │             │                      │
                             └────────────►│ SV-mutation          │
                                           │ enrichment (Fig 4) ⭐│
                                           │ temporal concordance │
                                           │ dose stratification  │
                                           │ (Fig 5, S3, S4)      │
                                           └──────────────────────┘
```

**The mutation half produces `*_dose_*_merged.csv` tables that the SV half
consumes** for breakpoint-proximal enrichment, temporal concordance, and
dose-stratified analyses. The mutation half is otherwise self-contained.

---

## Folder map

| Folder | What it covers | Manuscript sections |
|--------|----------------|---------------------|
| [`final codes_mutation/`](final%20codes_mutation/) | Point-mutation pipeline: VCF → SigProfiler → annotation → temporal patterns → Sankey → merged catalogues | §2.3, §2.5; Fig 2; Fig S1; Table S1 |
| [`final codes SV/`](final%20codes%20SV/) | Structural-variant pipeline: Manta/AnnotSV → PASS filter → temporal patterns → breakpoint-proximal enrichment → dose stratification → functional annotation | §2.4–§2.9; Figs 3, 4, 5; Figs S2, S3, S4; Tables S2, S3 |

Each folder contains its own `README.md`, `requirements.txt`, and (for the SV
folder) install notes for the external variant callers. Run the per-folder
READMEs in order; this top-level README only describes how the two halves
connect.

---

## Quick start

### 1. Install Python dependencies

The two folders have overlapping but not identical dependency lists. Install
them in one virtual environment:

```bash
python -m venv env
source env/bin/activate
pip install -r "final codes_mutation/requirements.txt"
pip install -r "final codes SV/requirements.txt"
```

(Both files pin compatible versions of pandas, numpy, scipy, matplotlib, and
seaborn. The mutation side adds SigProfilerExtractor + intervaltree; the SV
side adds dask + requests.)

### 2. Install upstream variant callers

These are external bioinformatics tools, not Python packages. Install them
in their own conda environments:

| Tool | Purpose | Instructions |
|------|---------|--------------|
| Mutect2 (GATK 4.5) | Somatic SNV/indel calling | manuscript §2.3 (no install script bundled) |
| Funcotator (GATK 4.5) | Mutation annotation | manuscript §2.3 |
| **Manta 1.6.0** | Somatic SV calling | [`final codes SV/install_manta.txt`](final%20codes%20SV/install_manta.txt) |
| **AnnotSV 3.3** | SV annotation | [`final codes SV/install_annotsv.txt`](final%20codes%20SV/install_annotsv.txt) |

The two `install_*.txt` files contain the exact commands used in the manuscript
including the critical `-annotationMode full` flag for AnnotSV.

### 3. External data

`final codes SV/11_fetch_annotations.py` requires the gnomAD v2.1.1 LoF
metrics file:

```
gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz
```

Download from <https://gnomad.broadinstitute.org/downloads>.

### 4. Run the mutation half first

```bash
cd "final codes_mutation"
# follow that folder's README — produces *_dose_*_merged.csv files
```

### 5. Run the SV half

```bash
cd "../final codes SV"
# follow that folder's README — several scripts read the merged CSVs from step 4
```

---

## Manuscript figure / table → script lookup

| Item | Producing script | Folder |
|------|------------------|--------|
| Fig 1 | Mermaid source `fig1.md` | `final codes SV/` |
| Fig 2 (DBS Sankey) | `sankey_visualization.py` | `final codes_mutation/` |
| Fig 3 (SV landscape) | `03_SV_landscape.py` | `final codes SV/` |
| Fig 4 (INV-DBS coupling) | `05_sv_mutation_correlation.py` + `07_sv_mut_vizualization.py` | `final codes SV/` |
| Fig 5 (dose response) | `10_dose_stratified.py` + `11_fetch_annotations.py` + `12_categorise_genes.py` + `13_dose_based_visualize.py` | `final codes SV/` |
| Fig S1 (SNV/MNS/InDel Sankeys) | `sankey_visualization.py` | `final codes_mutation/` |
| Fig S2 (repeats) | `04_repeat_analysis.py` | `final codes SV/` |
| Fig S3 (distance decay) | `06_sv_type_specific_decay.py` | `final codes SV/` |
| Fig S4 (temporal concordance) | `08_temporal_concordance.py` + `09_temporal_concordance_viz.py` | `final codes SV/` |
| Table S1 (pattern counts) | side-output of `mutation_pattern_assignment.py` | `final codes_mutation/` |
| Table S2 (Poisson enrichment stats) | side-output of `05_sv_mutation_correlation.py` | `final codes SV/` |
| Table S3 (high-constraint genes) | side-output of `12_categorise_genes.py` | `final codes SV/` |

The cross-references are also restated (with full method-section detail) in
each folder's own README.

---

## Conventions used throughout

- **Sample naming.** `d0` = unirradiated control. Doses `dA`–`dE` =
  0.001 / 0.01 / 0.1 / 1 / 2 mGy/hr. Timepoints `W1`, `W2`, `W3` are weekly
  post-exposure samples; `W0` is the pre-exposure baseline and is treated as
  "absent" in the leading position of every temporal pattern.
- **Pattern strings.** Four characters, position = `[W0][W1][W2][W3]`.
  Treatment present = `T`, control present = `C`, both = `B`, absent = `0`.
  The seven radiation-specific patterns are
  `0T00, 00T0, 000T, 0TT0, 0T0T, 00TT, 0TTT`.
- **AnnotSV `Annotation_mode`.** SV-counting scripts filter to
  `Annotation_mode == 'full'` so a single SV touching N genes is counted once,
  not N times.
- **Inversion deduplication.** AnnotSV can report the same INV with reversed
  coordinates; the relevant scripts collapse these using sorted `(start, end)`
  keys before counting.
- **Dose bins.** Low = `A, B, C` (0.001–0.1 mGy/hr); High = `D, E` (1–2 mGy/hr).

---

## Citation

If you use this pipeline, please cite the manuscript and the upstream tools:

- **Mutect2** — Benjamin D. *et al.* Calling somatic SNVs and indels with
  Mutect2. *bioRxiv* (2019). doi:10.1101/861054
- **SigProfilerMatrixGenerator** — Bergstrom E. N. *et al.* *BMC Genomics*
  20, 685 (2019). doi:10.1186/s12864-019-6041-2
- **SigProfilerExtractor** — Islam S. M. A. *et al.* *Cell Genomics* 2,
  100179 (2022). doi:10.1016/j.xgen.2022.100179
- **Manta** — Chen X. *et al.* *Bioinformatics* 32, 1220–1222 (2016).
  doi:10.1093/bioinformatics/btv710
- **AnnotSV** — Geoffroy V. *et al.* *Bioinformatics* 34, 3572–3574 (2018).
  doi:10.1093/bioinformatics/bty304

For the analysis framework itself, cite the LUCID manuscript (Chopra et al.).

---

## Repository scope

This repository contains only the **analytical pipeline** described in the
manuscript Methods. The following are intentionally out of scope:

- Cell culture and irradiation protocol — see Weinberg et al. (companion paper)
- Read alignment, QC, and variant calling — these are external upstream steps;
  Manta + AnnotSV install scripts are bundled, but Mutect2 / Funcotator /
  minimap2 / TrimGalore are documented in the manuscript Methods §2.2–§2.4
  rather than as install scripts here.
- Raw sequencing data — available via the BioProject accession in the
  manuscript Data Availability statement.
