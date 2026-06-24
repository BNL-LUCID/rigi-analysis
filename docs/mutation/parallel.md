# Parallel Helpers

Two outer-loop wrappers parallelize the slowest pipeline steps without
modifying the underlying logic. Each emits the same output layout as
the sequential script — drop-in replacements when sequential runs are
too slow.

## `parallel_annotate.py` (wraps step 4)

Chunk-level parallelism for `mutation_annotation.py`. Splits the input
`all_<TYPE>_mutations.pkl` into N row-balanced chunks, runs N annotation
subprocesses concurrently, and concatenates the per-chunk outputs into
a single combined `all_<TYPE>_annotated.pkl`.

**When to use**: SNV (tens of millions of rows). Reduces ~6 hours
sequential to ~30–40 minutes on a multi-core box. DBS / MNS / ID are
usually fast enough sequentially.

**Command:**
```bash
rigi-analysis-run parallel_annotate \
    --input          processed_data/all_SNV_mutations.pkl \
    --annotation-dir annotations \
    --build          hg38 \
    --output         annotated_mutations/SNV/all_SNV_annotated.pkl \
    --workers        24 \
    --no-csv         # optional: skip CSV write for huge tables
```

**Notes:**
- `--output` must be a `.pkl` file path, not a directory. The CSV is
  written alongside (filename derived by replacing `.pkl` → `.csv`)
  unless `--no-csv` is passed.
- Each worker loads its own copy of the interval-tree annotation data
  (a few hundred MB to ~1 GB on hg38). RAM scales linearly with
  `--workers`.
- One failed chunk fails the whole run with exit 1. Per-chunk logs
  land in `<output>_parallel_work/out/chunk_NNN/subprocess.log`. The
  scratch directory is preserved on failure (recover via `pd.concat`
  of the surviving chunk pickles); deleted on success unless
  `--keep-work-dir` is passed.

## `parallel_pattern_assignment.py` (wraps step 5)

Dose-level parallelism for `mutation_pattern_assignment.py`. Splits
the annotated input by dose (each split contains the controls plus
one dose), launches one subprocess per non-control dose, each writing
to its own `dose_<dose>/` subdirectory. No merge step — each subprocess
produces its own final per-dose output.

**When to use**: any mutation type where sequential pattern assignment
is too slow. Capped at the number of non-control doses in the data
(typically 5).

**Command:**
```bash
rigi-analysis-run parallel_pattern_assignment \
    --input          annotated_mutations/SNV/all_SNV_annotated.pkl \
    --output-dir     pattern_analysis/SNV \
    --mutation-type  SNV \
    --workers        5
```

**Notes:**
- Output layout is identical to a sequential run. Verify on a small
  type first by diffing against an existing sequential output:
  ```bash
  diff <(sort old/dose_dA/mutation_annotations_dose_dA_DBS.csv) \
       <(sort new/dose_dA/mutation_annotations_dose_dA_DBS.csv)
  # expected: no output
  ```
- Memory caveat: each worker loads its filtered DataFrame (controls
  plus one dose ≈ 30% of the full input). 5 workers × 30% ≈ 1.5×
  peak memory vs sequential. On large inputs (SNV) this can OOM-kill
  workers — kernel sends SIGKILL, workers die with returncode -9 and
  empty logs (no chance to flush). If you see this pattern, drop
  `--workers` to 2 or 3, or run failed doses sequentially using the
  per-dose split inputs in the work directory (preserved on failure).
- The script captures each subprocess's output into
  `<output-dir>_parallel_work/logs/<dose>.log`. Empty logs after a
  SIGKILL is the OOM-kill signature.

## Caveats common to both helpers

- The wrappers don't modify the underlying scripts, so all bug fixes
  and guards in `mutation_annotation.py` and
  `mutation_pattern_assignment.py` apply inside each subprocess.
- Per-chunk QC plots from `mutation_annotation.py` (the
  `quality_transcription_analysis/` and `indel_mechanism_analysis/`
  directories) are produced per chunk and discarded with the scratch
  directory — they're not representative when computed on a slice. If
  you need global QC plots after parallel annotation, the merged
  DataFrame is what matters; the QC plots are diagnostic only.
