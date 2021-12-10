#!/usr/bin/env zsh
#SBATCH --gres=gpu:1
#SBATCH --time=3-0:0:0
#SBATCH --partition=research
#SBATCH -o OUTPUT/tumbler-%A_%a.out
#SBATCH -e OUTPUT/tumber-%A_%a.err
##SBATCH --account=sbel
##SBATCH --qos=sbel_owner
##SBATCH --qos=priority

##SBATCH --array=1-16
module load nvidia/cuda/11.3.1

drum_height=3
drum_omega=10

cd ../build
./test_GPU_tumbler_settling ${drum_height}
./test_GPU_tumbler_spinning ${drum_height} ${drum_omega}


