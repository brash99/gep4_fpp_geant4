#!/bin/bash
#SBATCH --job-name=gep4_fpp_6p7
#SBATCH --account=halla
#SBATCH --partition=production
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/farm_out/brash/gep4_fpp_mono_6p7_%j.out
#SBATCH --error=/farm_out/brash/gep4_fpp_mono_6p7_%j.err

set -euo pipefail

source /etc/profile.d/modules.sh
module purge
module use /group/halla/modulefiles
module load geant4/11.2.1
module load root/6.30.04

repo=/work/hallc/gep/brash/gep4_fpp_geant4
run_dir="$repo/runs/mono_6p7"
mkdir -p "$repo/logs" "$run_dir"

cd "$repo"
python3 scripts/generate_mono_sample.py \
  "$run_dir/events.csv" --events 100000 --momentum 6.7

./build/gep4_fpp_sim \
  macros/run_100k.mac \
  "$run_dir/events.csv" \
  "$run_dir/transport.root" \
  FTFP_BERT

python3 analysis/summarize_efficiency.py \
  "$run_dir/transport.root" --theta-min 5 --theta-max 30 \
  > "$run_dir/efficiency_5_30deg.txt"
