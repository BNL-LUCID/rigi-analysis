# RIGI Analysis: Spatio-Temporal Analysis of Radiation-Induced Genomic Instability

Analytical code accompanying:

> **A Spatio-Temporal Analysis Framework for Characterizing 
> Radiation-Induced Genomic Instability** (2026) bioRxiv 
> [2026.02.21.707188](https://doi.org/10.64898/2026.02.21.707188)

The pipeline analyses whole-genome sequencing of HUVEC cells exposed 
to chronic low-dose gamma radiation (0.20–2.62 mGy/hr) over three 
weeks. It produces the temporal pattern catalogues, breakpoint-proximal 
mutation enrichment statistics, dose-stratified analyses, and figures 
presented in the paper.

## Repository layout

```
src/rigi_analysis/
├── mutation/       # Point-mutation pipeline (Mutect2 → SigProfiler → patterns → Sankey)
├── sv/             # Structural-variant pipeline (Manta → AnnotSV → enrichment → figures)
├── workflows/      # Async workflow orchestrators (MutationPipeline, SVPipeline, FullWorkflow)
└── utils/          # Shared utilities (CLI dispatcher, datetime helpers)
docs/
├── mutation/       # Mutation pipeline documentation
└── sv/             # SV pipeline documentation
tests/              # Dry-run workflow tests (no HPC required)
```

The two pipelines meet at an intermediate artefact: dose-merged mutation
tables (`*_dose_*_merged.csv`) produced by the `mutation` half and consumed
by the `sv` half. **Run the `mutation` pipeline first, then the `sv` pipeline, 
or run the full workflow only, which will do both for you.**

## Setup

```bash
conda env create -f environment.yml
conda activate rigi-analysis
pip install -e ".[dev]"
```

For GPU-accelerated SigProfiler runs, install the optional `gpu` extra:

```bash
pip install -e ".[gpu]"
```

Full documentation, including external tool setup (`Mutect2`, `Manta`, 
`AnnotSV`) and pipeline run instructions, is in [`docs/`](docs/README.md).

## Citation

If you use this pipeline, please cite the manuscript above. Upstream 
tools (`Mutect2`, `SigProfiler`, `Manta`, `AnnotSV`) should also be cited 
separately — see the manuscript References section.

## Funding

This work was supported by the U.S. Department of Energy, Office of 
Science, Biological and Environmental Research program 
(B&R# KP1601017, FWP# CC140).

## Out of scope

Raw sequencing data, alignment, and upstream variant calling are not 
included. See the manuscript Methods (§4) and the BioProject accession 
in Data Access for those.
