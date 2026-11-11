#!/bin/bash
#SBATCH --job-name=gep4_q4_layout
#SBATCH --partition=production
#SBATCH --account=halla
#SBATCH --array=0-3%2
#SBATCH --time=03:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/farm_out/brash/gep4_q4_layout_%A_%a.out
#SBATCH --error=/farm_out/brash/gep4_q4_layout_%A_%a.err

set -euo pipefail
source /etc/profile.d/modules.sh
module purge
module use /group/halla/modulefiles
module load geant4/11.2.1
module load root/6.30.04

repo=/work/hallc/gep/brash/gep4_fpp_geant4
tags=(legacy_55p88 compact_55p88 compact_85 compact_90)
thicknesses=(55.88 55.88 85 90)
tag=${tags[$SLURM_ARRAY_TASK_ID]}
thickness=${thicknesses[$SLURM_ARRAY_TASK_ID]}
run_dir="$repo/runs/q4_${SLURM_ARRAY_JOB_ID}/${tag}"
mkdir -p "$run_dir"
cd "$repo"

python3 scripts/generate_mono_sample.py "$run_dir/events.csv" \
  --events 100000 --momentum 6.7 --seed 20261105

layout_args=()
if [[ "$tag" == compact_* ]]; then
  layout_args+=(--compact-downstream-fpp)
fi

./build-q4/gep4_fpp_sim macros/run_100k.mac "$run_dir/events.csv" \
  "$run_dir/transport.root" --physics FTFP_BERT --single-analyzer \
  --analyzer-thickness-cm "$thickness" "${layout_args[@]}"

python3 analysis/compare_track_selection.py "$run_dir/transport.root" \
  "$run_dir/track_selection" --thickness "$thickness" \
  > "$run_dir/track_selection.txt"
python3 analysis/plot_track_selection.py "$run_dir/track_selection"
python3 analysis/calculate_fom.py "$run_dir/transport.root" \
  "$run_dir/fom" --thickness "$thickness" \
  > "$run_dir/fom.txt"
python3 analysis/analyze_four_chamber_reconstruction_closure.py \
  "$run_dir/transport.root" "$run_dir/all4_closure" \
  --thickness "$thickness" > "$run_dir/all4_closure.txt"

for min_angle in 0.5 1.0 1.5 2.0; do
  angle_tag=${min_angle/./p}
  python3 analysis/compare_track_selection.py "$run_dir/transport.root" \
    "$run_dir/track_selection_min_${angle_tag}" --thickness "$thickness" \
    --min-angle-deg "$min_angle" \
    > "$run_dir/track_selection_min_${angle_tag}.txt"
done

echo "$tag complete"
