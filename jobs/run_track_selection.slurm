#!/bin/bash
#SBATCH --job-name=fpp_tracks
#SBATCH --partition=production
#SBATCH --account=halla
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/farm_out/brash/fpp_tracks_%j.out
#SBATCH --error=/farm_out/brash/fpp_tracks_%j.err
set -euo pipefail
source /etc/profile.d/modules.sh
module use /group/halla/modulefiles
module load geant4/11.2.1 root/6.30.04
repo=/work/hallc/gep/brash/gep4_fpp_geant4
cd "$repo"
run_dir="$repo/runs/track_selection_${SLURM_JOB_ID}"
mkdir "$run_dir"
git rev-parse HEAD > "$run_dir/source_head.txt"
git diff > "$run_dir/source.patch"
sha256sum build-track-selection/gep4_fpp_sim runs/single_55p88cm/events.csv > "$run_dir/input_hashes.txt"
./build-track-selection/gep4_fpp_sim macros/run_100k.mac \
  runs/single_55p88cm/events.csv "$run_dir/transport.root" \
  --physics FTFP_BERT --single-analyzer --analyzer-thickness-cm 55.88 \
  > "$run_dir/transport.log" 2>&1
python3 analysis/compare_track_selection.py "$run_dir/transport.root" \
  "$run_dir/comparison" --baseline-root runs/single_55p88cm/transport.root \
  > "$run_dir/comparison.log" 2>&1
python3 analysis/plot_track_selection.py "$run_dir/comparison"
tar -czf "$run_dir/source_snapshot.tar.gz" CMakeLists.txt include/EventAction.hh \
  src/EventAction.cc src/RunAction.cc analysis/track_selection.py \
  analysis/compare_track_selection.py analysis/plot_track_selection.py \
  analysis/analyzing_power.py analysis/calculate_fom.py \
  tests/test_track_selection.py jobs/run_track_selection.slurm TRACK_SELECTION_STUDY.md
