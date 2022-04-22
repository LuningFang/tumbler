#!/usr/bin/env zsh
#SBATCH --gres=gpu:1
#SBATCH --time=3-0:0:0
#SBATCH --partition=research
#SBATCH -o nozzle-%A.out
#SBATCH -e nozzle-%A.err
##SBATCH --account=sbel
##SBATCH --qos=sbel_owner
##SBATCH --qos=priority

##SBATCH --array=1-16
module load nvidia/cuda/11.3.1

cd ../build

./test_GPU_nozzle \
--nozzle_diameter="0.1" \
--nozzle_angle="60" \
--output_directory="nozzle" \
--mu_r="0.1" \
--mu_s="0.15"

