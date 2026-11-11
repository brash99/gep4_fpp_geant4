# Single-track and smallest-angle FPP comparison

This study addresses Andrew Puckett's November 10, 2026 single-track question.
The published 26.264% useful-pT efficiency selected a leading-momentum,
truth-associated charged track, with no event-level multiplicity veto.

## Definitions

The new `candidates` ROOT tree retains one row per GEANT4 track ID crossing
both FPP1 and FPP2, before any cone, closest-approach, or angular cuts. Neutral
particles are excluded; charged primary and secondary particles are included.
Repeated crossings of one track retain the existing crossing-map behavior.
A candidate's chamber mask identifies whether it also crosses FPP3 and FPP4.
Multiplicity is a geometric proxy for reconstructibility, not measured
pattern-recognition multiplicity.

Four policies use the same transport sample and incident-proton denominator:

- `legacy`: preserve the original independent highest-momentum candidate
  choices for the four- and two-chamber topologies, solely as a reference.
- `single`: require exactly one front-pair candidate before any analysis cuts.
  That candidate must cross all four chambers below 12 degrees. This also
  vetoes extra front-pair tracks that miss the downstream pair.
- `min_angle_uncut`: retain the previous smallest-angle rule with no explicit
  one-degree lower bound, as a regression reference.
- `min_angle`: first restrict front-pair tracks to polar angle strictly greater
  than 1 degree, then choose the smallest angle, with track ID as a deterministic
  tie breaker. The bound is configurable with `--min-angle-deg` (default 1).
  A track exactly at 1 degree is excluded. Tracks below the bound do not prevent
  another track above it from being chosen. Require the chosen track to cross
  all four chambers below 12 degrees. No outgoing momentum is used in the
  choice. A track failing a remaining cut is not replaced by another candidate.

After choosing, apply the existing cone test, analyzer z-close interval,
s-close <= 1 cm, and the 0.07 <= pT <= 1.0 GeV/c window. The analyzer remains
55.88 cm, incident momentum 6.7 GeV/c, and physics list FTFP_BERT. Incoming
proton momentum, known from the upstream spectrometer in the intended
experiment, enters pT and Ay; it is distinct from the unmeasured outgoing
track momentum that was used in the legacy ranking.

The minimum-angle policy includes every accepted single-track event whose
track is strictly above its angular threshold. It may
lose events relative to the legacy policy when its chosen track fails cuts or
lies below the useful-pT region. Per-event choices are saved for auditing.

## Interpretation limits

Angles and closest-approach quantities remain ideal truth-associated values,
as in the original FOM. The existing unsigned-drift wire fits are preserved
as closure diagnostics; these selections do not measure real multi-track
pattern recognition or chamber inefficiencies. The empirical GEp-III Ay
model is a single-track calibration. Applying it to recovered multi-track
events is an explicitly provisional FOM estimate, not a validation of that
sample's analyzing power or statistical independence.

## Reproduction and validation

Build in a separate directory to preserve the existing executable:

```sh
cmake -S . -B build-track-selection -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DPython3_EXECUTABLE="$(command -v python3)"
cmake --build build-track-selection -j4
ctest --test-dir build-track-selection --output-on-failure
sbatch jobs/run_track_selection.slurm
```

The job reuses the original 100,000 incident-proton CSV and does not overwrite
old ROOT files. It records HEAD, the working diff, and input/binary hashes.
The analysis records source hashes, checks candidate uniqueness and source
event IDs, reconstructs the original selection from the candidate table,
compares all original physics branches against the old tree by entry, and
checks that the single-track useful sample above the angular threshold is
contained in the minimum-angle sample. JSON and summary CSV files contain efficiencies, RMS effective Ay,
intrinsic and beam-polarization-weighted FOM, and Monte Carlo standard errors.
The audit CSV records one row per incident event including zero-track events.

The old 100,000-event ROOT file has event_id=-1, weight=1, incident_p_GeV=0
because it predates the input-metadata persistence correction already present
in HEAD. The analysis explicitly checks these defaults in every old entry
and excludes those three metadata fields from the physics-branch comparison.
The new run retains its actual CSV event IDs, weights, and 6.7 GeV/c momentum.
The old published FOM used the monoenergetic momentum override.

Use the ROOT-compatible Python selected by the farm environment (Python 3.9
for ROOT 6.30.04). CMake otherwise discovers Python 3.13, whose ROOT
bindings are unavailable. This was diagnosed by the integration check and
resolved through the explicit interpreter configuration above.

## Completed 100,000-event comparison

Job 14574698 completed normally on November 10, 2026. The original physics
branches agree exactly with the previous sample after the documented metadata
exclusions, and the original selection is reproduced from 240,057 candidate
rows. Outputs are in `runs/track_selection_14574698/`.

| Policy | Useful efficiency | Polarized FOM / incident, Pe=0.85 |
|---|---:|---:|
| Original leading momentum | 26.264% | 0.000755832045 |
| Exclusive single front-pair track | 12.492% | 0.000309876813 |
| Previous smallest angle, no 1-degree bound | 28.939% | 0.000851083403 |

The smallest-angle sample adds 16,447 useful multi-track events to the 12,492
single-track events. Its provisional FOM is 12.60% above the original result,
but the single-track Ay calibration has not been validated for that recovered
population. See the run's `REPORT.md`, JSON/CSV summaries, per-event audit,
plots, and source snapshot for definitions, checks, and interpretation limits.

## Updated default: smallest angle strictly above 1 degree

The default now filters theta <= 1 degree before ranking. Reanalysis of the
same saved transport produces 25,892 useful events (25.892%), effective
Ay=0.06704, intrinsic FOM=0.00116370570 per incident proton, and polarized
FOM=0.000840777365 at Pe=0.85. This retains 98.79% of the previous
smallest-angle FOM. The new sample comprises 9,663 useful single-track events
above one degree and 16,229 multi-track events. The unchanged single-track
reference still has no explicit one-degree cut. Updated outputs use the
comparison_gt1 stem and REPORT_gt1.md; earlier outputs are preserved.
