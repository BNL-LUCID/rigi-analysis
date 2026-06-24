# Pattern Encoding Reference

Compact notation used throughout the mutation pipeline to describe a
mutation's presence/absence trajectory across the four sequencing
timepoints (W0–W3) and two conditions (control, treated).

## Categorical states

Per timepoint, a mutation is in one of four states based on which arms it
appears in:

| Symbol | Meaning |
|--------|---------|
| `C`    | Detected in control samples only |
| `T`    | Detected in treated (irradiated) samples only |
| `B`    | Detected in both control and treated samples |
| `0`    | Absent from both arms |

A four-position string `[W0][W1][W2][W3]` encodes the full trajectory.
Example: `0T0T` means absent at W0, treated-only at W1, absent at W2,
treated-only at W3 (a recurrent radiation-induced mutation).

## Example patterns

| Pattern | Interpretation |
|---------|----------------|
| `CBBB`  | Present at baseline, persists in both control and treated arms throughout |
| `0T00`  | Treatment-specific, W1 only |
| `0TTT`  | Treatment-induced, persists across all post-exposure timepoints |
| `00TT`  | Treatment-induced, appears at W2 and persists to W3 |
| `0T0T`  | Treatment-specific, recurrent (W1, absent W2, reappears W3) |
| `C000`  | Control-specific at baseline, lost thereafter |

## Pattern categories

Patterns are grouped into four high-level categories used by downstream
scripts (`Pattern_Group` column in merged outputs):

- **Radiation-specific** — `0T00`, `00T0`, `000T`, `0TT0`, `00TT`, `0T0T`, `0TTT`
- **Control-specific** — `C000`, `0C00`, `00C0`, `000C`
- **Baseline** — `CBBB` (present at all timepoints in both conditions)
- **Other** — any other combination

## Sankey state taxonomy

The Sankey scripts (`compute_sankey.py` / `render_sankey.py`) use an
expanded state taxonomy that distinguishes new appearances from
reappearances:

- **W0**: `Present`, `Absent`
- **W1**: `Both`, `Lost`, `Control`, `Exposed` (no Recurrent — W0 is baseline)
- **W2**: `Both`, `Exposed_Recurrent`, `Control_Recurrent`, `Exposed`, `Control`, `Lost`
- **W3**: same as W2

A state is `*_Recurrent` only if the mutation was present in that arm at
some prior timepoint AND was absent in that arm at the immediately
preceding timepoint. Newly appearing mutations without prior presence
are classified as `Exposed`/`Control`, not Recurrent. This taxonomy is
what the manuscript Fig 2 visualization uses.

## Dose and timepoint conventions

- **Doses**: `Control` (unirradiated), `dA`–`dE` (radiation, ascending dose rate)
- **Timepoints**: `W0` (pre-exposure baseline), `W1`, `W2`, `W3` (weekly post-exposure)
- `W0` is treated as baseline-absent in trajectory leading position regardless of
  whether the mutation was detected in the W0 control sample.
