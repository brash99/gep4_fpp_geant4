# GEp-4 Hall C FPP GEANT4 simulation

Standalone GEANT4 transport simulation for estimating the Hall C focal-plane
polarimeter efficiency and figure of merit for a possible GEp-4 reboot.

The initial geometry is transcribed from the GEp-III GEANT3 program at:

```text
/work/hallc/gep/brash/e93027/brash/geant_gep
```

This repository intentionally remains separate from SIMC. SIMC supplies
weighted, HMS-accepted proton events at the FPP reference plane; this program
transports them through the downstream scintillators, CH2 analyzers, and
idealized FPP tracking planes.

## Current scope

- GEp-III baseline geometry with two 55.88 cm CH2 analyzers.
- CSV input for weighted proton phase space.
- `FTFP_BERT` hadronic physics by default.
- Idealized FPP tracking planes.
- ROOT ntuple containing incident kinematics, first hadronic interaction,
  scattering angle, outgoing momentum, and per-track chamber-pair acceptance.

The current output is a transport-validation product, not yet a final FOM.
Analyzing power will be applied outside GEANT4 after defining and validating
the useful-scattering selection.

## Input contract

CSV columns are:

```text
event_id,weight,x_cm,y_cm,z_cm,px_GeV,py_GeV,pz_GeV
```

Coordinates refer to the simulation reference frame. The GEp-III baseline
uses the HMS focal plane as `z = 0`; downstream elements retain the positions
encoded in the legacy program.

## Build and run on ifarm

```bash
module load geant4/11.2.1
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j4
ctest --test-dir build --output-on-failure
```

Run explicitly with:

```bash
./build/gep4_fpp_sim macros/run_example.mac data/example_events.csv output.root
```

An optional fourth argument selects a Geant4 reference physics list, for
example `QGSP_BERT`.

Summarize the weighted geometry-level efficiency for a chosen useful angular
range with:

```bash
python3 analysis/summarize_efficiency.py output.root --theta-min 5 --theta-max 30
```

The reported statistical uncertainty is only a weighted-sample diagnostic.
Hadronic-model, geometry, reconstruction, and analyzing-power uncertainties
must be evaluated separately.

The primary reconstruction diagnostic is the two-dimensional distribution of
polar scattering angle versus the z coordinate of closest approach of the
upstream and downstream tracks:

```bash
python3 analysis/plot_theta_vs_zclose.py output.root theta_vs_zclose
```

The plot marks the physical z boundaries of each CH2 analyzer. This is the
direct analogue of the central GEp-III FPP reconstruction diagnostic.
The same command also writes a companion canvas containing
`theta_polar` versus `s_close` and the accepted-event azimuthal distributions
for both analyzers. The azimuth follows the legacy GEp-III convention:
`atan2(delta horizontal slope, delta vertical slope)`, stored in radians and
wrapped onto `[0, 2 pi)`.

Scan the geometry-level efficiency after requiring `z_close` to lie inside
the corresponding analyzer and varying the `s_close` requirement with:

```bash
python3 analysis/scan_reconstruction_cuts.py output.root \
  --theta-min 5 --theta-max 30
```

The combined efficiency is an event-level union: an event accepted in both
analyzers is counted only once.

## Monoenergetic validation sample

Generate a deterministic central-ray sample with:

```bash
python3 scripts/generate_mono_sample.py events.csv \
  --events 100000 --momentum 6.7
```

Optional Gaussian position and projected-angle widths support controlled
acceptance scans. The committed Slurm job `jobs/run_mono_6p7.slurm` runs the
first 100,000-event central-ray benchmark used for GEANT3/GEANT4 closure.

## Validation plan

1. Close the geometry against the GEp-III source and survey.
2. Compare monoenergetic proton transport between GEANT3 and GEANT4.
3. Reproduce measured GEp-III/GEp-2gamma efficiencies and angular spectra.
4. Process weighted GEp-4 SIMC events.
5. Apply interchangeable analyzing-power parameterizations externally.

See the parent study's `fpp_geant3_to_geant4_assessment.md` for the detailed
migration rationale.
