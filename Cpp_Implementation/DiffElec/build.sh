# ---------------------------------------- #
#         DiffElec build script            #
# ---------------------------------------- #

# ----------------------------------------- #
# This script simply calls CMake from the   #
# build directory and compiles the library. #
#                                           #
# Usage: bash build.sh {make_arg}           #
# ----------------------------------------- #

# Change to file directory. 
cd "$(dirname "$(realpath "$0")")";

# Check if build/ dir exists.
if [ ! -d build ]; then
    mkdir build
else
    rm -rf build
    mkdir  build
fi

# Change to build dir and compile the library.
cd build
cmake ..
make $@
