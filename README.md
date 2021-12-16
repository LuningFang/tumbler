# tumbler
source code and scripts for Metal Powder Group to run tumbler tests

## Things to have ready on your own machine
WinSCP for Windows user to download/edit files from the remote location

Paraview for visualize particles

## Get a copy of Chrono and its submodules
0. ssh into your account on Euler
````sh
ssh [account_name]@euler-login-2.wacc.wisc.edu
````
1. Clone chrono repo
````sh
git clone https://github.com/projectchrono/chrono.git
````
2. Now you have a copy of chrono, go into chrono directory
````sh
cd chrono
````
3. set up submodules by typing the following commands
````sh
git submodule init
git submodule update
````

## Build chrono
1. go to the home directory, make a new folder and go into the new folder
````sh
cd
mkdir build_chrono
cd build_chrono
````
Here are some 
2. load modules needed for building Chrono::GPU
````sh
module load nvidia/cuda/11.3.1.lua 
````
3. load cmake, first argument locates the source code, second argument is the build directory (in this case current folder .), third argument is a compiler flag
````sh
ccmake ../chrono . -DCUDA_HOST_COMPILER=$CU_CC
````
4. In cmake gui, type `c` once for configure, scroll down to `ENABLE_MODULE_GPU` and set to `ON`, press `c` twice, and then `g` for generate to create `Makefile`. If there is no error, cmake gui will close, are you are back to `build_chrono` directory. Now type
````
make -j 64
````
it will take its sweet time, so be patient. Note that you do not need to modify anything in `chrono` or `build_chrono`  

## Get a copy of tumber test repo and compile the code 
1. go to your home directory by typing `cd`, then clone tumbler test repo
````sh
git clone https://github.com/LuningFang/tumbler.git
````
Folder `source` has two programs, `test_GPU_tumbler_settling.cpp` is used to create a settled configuration of particles and `test_GPU_tumbler_spinning.cpp` is for rotating the tumbler. Folder `scripts` contains the script necessary for running the above two tests  
2. build the tests, create the build folder
````sh
cd tumbler
mkdir build
cd build
````
use cmake to generate the `Makefile`
````sh
ccmake ../source/ .
````
you need to link tumbler test to Chrono library, type `c` once for configure, then type `e` to exit screen, then in the field of `Chrono_DIR` type the following directory `srv/home/[YOUR_EULR_ACCOUNT_NAME]/build_chrono/cmake/`, type `c` twice until `g` shows up, press `g`  
Now you can build
````sh
make -j 64
````
## Submit jobs to run tumbler tests
1. Go to the directory of scripts
````sh
cd ../scripts/
````
2. Submit `run_tumbler_test.sh` using the following command
````sh
sbatch run_tumbler_test.sh
````
To see where your job is and how long it has been running, you can type the following command in the terminal, 
````sh
squeue -u [YOUR_EULER_ACCOUNT_NAME]
````
3. In `run_tumbler_test.sh` you can modify parameter `drum_height` and `drum_omega`, note that drum height is in centimeter and drum omega is in rpm. If your job runs successfully, you will see an output file named as `tumbler*.out` in the same folder. Each time step, a csv file that contains the positions and absolute velocity of all the particles is created and written in the following folder 
`/srv/home/fang/tumbler/build/DEMO_OUTPUT/tumbler_settling/` for settling phase, and   `/srv/home/fang/tumbler/build/DEMO_OUTPUT/tumbler_spinning/spinning_omega_XX_rpm/` for spinning phase at XX rpm. 
4. To modify `drum_height` and `drum_omega` in the script, Windows users can use WinSCP, Mac users can use `vim`, see [this link](https://www.youtube.com/watch?v=ggSyF1SVFr4) for a tutorial video.

## Postprocessing
1. After the simulation is done, you can download your output from Euler to local machine. 

Windows user can use WinSCP.   
Hostname is euler-login-2.wacc.wisc.edu. 
Username and password is your username and password for Euler.   

For Mac user, in your terminal, go to any directory where you wish to keep the files, type
````sh
scp '[YOUR_EULER_ACCOUNT_NAME]@euler-login-2.wacc.wisc.edu:[directory-to-your-csv-files]/*' .
````
2. Start Paraview for postprocessing, slides available

## What to do when Luning modifies the code
1. Go to the tumbler directory on your cluster using `cd`, use command `pwd` to check if you are in the right place, my output look something like this `srv/home/fang/tumbler`
2. Type the command to get updated source code
````sh
git pull
````
3. Rebuild the source code to get an updated executable. Go to the build directory (use `pwd` to check, mine looks like `/srv/home/fang/tumbler/build`), and type
````sh
make -j 64
````
Let Luning know if any of the above two steps don't work  
4. Switch to your script folder,
````sh
cd ../scripts
````
and run the script using `sbatch` command. 
