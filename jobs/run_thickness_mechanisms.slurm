#!/bin/bash
#SBATCH --job-name=gep4_thickness_mechanisms
#SBATCH --partition=production
#SBATCH --account=halla
#SBATCH --array=0-4
#SBATCH --time=00:45:00
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1
#SBATCH --output=/farm_out/brash/gep4_thickness_mechanisms_%A_%a.out
#SBATCH --error=/farm_out/brash/gep4_thickness_mechanisms_%A_%a.err

set -euo pipefail
source /etc/profile.d/modules.sh
module purge
module use /group/halla/modulefiles
module load geant4/11.2.1
module load root/6.30.04

repo=/work/hallc/gep/brash/gep4_fpp_geant4
thicknesses=(30 55.88 70 85 90)
tags=(compact_30 compact_55p88 compact_70 compact_85 compact_90)
thickness=${thicknesses[$SLURM_ARRAY_TASK_ID]}
tag=${tags[$SLURM_ARRAY_TASK_ID]}
run_dir="$repo/runs/thickness_mechanisms_${SLURM_ARRAY_JOB_ID}/${tag}"
mkdir -p "$run_dir"
cd "$repo"

python3 scripts/generate_mono_sample.py "$run_dir/events.csv" \
  --events 100000 --momentum 6.7 --seed 20261111
./build-mechanisms/gep4_fpp_sim macros/run_100k_seeded.mac \
  "$run_dir/events.csv" "$run_dir/transport.root" \
  --physics FTFP_BERT --single-analyzer --compact-downstream-fpp \
  --analyzer-thickness-cm "$thickness"

python3 analysis/compare_track_selection.py "$run_dir/transport.root" \
  "$run_dir/track_selection" --thickness "$thickness" \
  > "$run_dir/track_selection.txt"
python3 analysis/analyze_thickness_mechanisms.py \
  --sample "$tag" "$thickness" "$run_dir/transport.root" \
  --output-stem "$run_dir/mechanisms" > "$run_dir/mechanisms.txt"
python3 analysis/analyze_four_chamber_reconstruction_closure.py \
  "$run_dir/transport.root" "$run_dir/all4_closure" \
  --thickness "$thickness" > "$run_dir/all4_closure.txt"

echo "$tag complete"
