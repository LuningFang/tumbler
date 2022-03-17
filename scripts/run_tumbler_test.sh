#!/usr/bin/env zsh
#SBATCH --gres=gpu:1
#SBATCH --time=3-0:0:0
#SBATCH --partition=research
#SBATCH -o tumbler-%A.out
#SBATCH -e tumbler-%A.err
##SBATCH --account=sbel
##SBATCH --qos=sbel_owner
##SBATCH --qos=priority

##SBATCH --array=1-16
module load nvidia/cuda/11.3.1

## drum height in cm, omega in rpm
drum_height=0.3
drum_omega=15

cd ../build

./test_GPU_tumbler_settling \
--particle_radius="0.05" \
--output_directory="some_test" \
--mu_r="0.1" \
--mu_s="0.3" \
--step_size="0.000005" \
--drum_height=${drum_height}

 ./test_GPU_tumbler_spinning \
 --particle_radius="0.05" \
 --output_directory="some_test" \
 --mu_r="0.1" \
 --mu_s="0.3" \
 --step_size="0.000005" \
 --drum_height="0.3" \
 --drum_speed=${drum_omega}
