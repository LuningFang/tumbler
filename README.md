# tumbler
source code and scripts for Metal Powder Group to run tumbler tests

## Get a copy of Chrono and its submodules
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
you need to link tumbler test to Chrono library, type `c` once for configure, in the field of `Chrono_DIR` type the following directory `srv/home/[YOUR_EULR_ACCOUNT_NAME]/chrono_build/cmake/`, type `c` twice until `g` shows up, press `g`  
Now you can build
````sh
make -j 64
````
