# LUCID — Spatio-Temporal Analysis of Radiation-Induced Genomic Instability

Analytical code accompanying:

> **A Spatio-Temporal Analysis Framework for Characterizing 
> Radiation-Induced Genomic Instability**
> Chopra K., Cucinell C., Weinberg R., Forrester S., Brettin T., 
> Kilic O. O., Yoon B. J. (2026).

The pipeline analyses whole-genome sequencing of HUVEC cells exposed 
to chronic low-dose gamma radiation (0.20–2.62 mGy/hr) over three 
weeks. It produces the temporal pattern catalogues, breakpoint-proximal 
mutation enrichment statistics, dose-stratified analyses, and figures 
presented in the paper.

## Repository layout

The pipeline has two halves that meet at one intermediate artefact:

- **`final codes_mutation/`** — point-mutation pipeline (Mutect2 → 
  SigProfiler → temporal patterns → Sankey). Produces 
  `*_dose_*_merged.csv` tables.
- **`final codes SV/`** — structural-variant pipeline (Manta → AnnotSV 
  → temporal patterns → breakpoint enrichment → dose stratification). 
  Reads the merged CSVs from the mutation half.

Run the mutation half first, then the SV half. Each folder has its own 
README with run instructions, figure/table cross-references, and 
dependencies.

## Setup

```bash
python -m venv env
source env/bin/activate
pip install -r "final codes_mutation/requirements.txt"
pip install -r "final codes SV/requirements.txt"
```

External tools (Mutect2, Manta, AnnotSV) are documented in the per-folder 
READMEs and `install_*.txt` notes.

## Citation

If you use this pipeline, please cite the manuscript above. Upstream 
tools (Mutect2, SigProfiler, Manta, AnnotSV) should also be cited 
separately — see the manuscript References section.

## Funding

This work was supported by the U.S. Department of Energy, Office of 
Science, Biological and Environmental Research program 
(B&R# KP1601017, FWP# CC140).

## Out of scope

Raw sequencing data, alignment, and upstream variant calling are not 
included. See the manuscript Methods (§4) and the BioProject accession 
in Data Access for those.
