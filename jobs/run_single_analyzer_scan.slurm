#!/bin/bash
#SBATCH --job-name=gep4_fpp_single
#SBATCH --partition=production
#SBATCH --account=halla
#SBATCH --array=0-3%2
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/farm_out/brash/gep4_fpp_single_%A_%a.out
#SBATCH --error=/farm_out/brash/gep4_fpp_single_%A_%a.err

set -euo pipefail
source /etc/profile.d/modules.sh
module purge
module use /group/halla/modulefiles
module load geant4/11.2.1
module load root/6.30.04

repo=/work/hallc/gep/brash/gep4_fpp_geant4
thicknesses=(40 50 55.88 60)
thickness=${thicknesses[$SLURM_ARRAY_TASK_ID]}
tag=${thickness/./p}
run_dir="$repo/runs/single_${tag}cm"
mkdir -p "$run_dir"
cd "$repo"

python3 scripts/generate_mono_sample.py "$run_dir/events.csv" \
  --events 100000 --momentum 6.7
./build/gep4_fpp_sim macros/run_100k.mac "$run_dir/events.csv" \
  "$run_dir/transport.root" --physics FTFP_BERT --single-analyzer \
  --analyzer-thickness-cm "$thickness"
python3 analysis/summarize_single_analyzer.py "$run_dir/transport.root" \
  --thickness "$thickness" --theta-min 5 --theta-max 30 --sclose-max 1 \
  > "$run_dir/efficiency.txt"
