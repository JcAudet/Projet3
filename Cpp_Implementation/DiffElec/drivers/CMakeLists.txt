# ------------------------------------------------------------------------------------------- #
# Author:       Jean-Christophe Gagnon-Audet  <Jean-Christophe.gagnon-Audet@emt.inrs.ca>      #
# Date:         2019 - 05 - 15                                                                #
# Description:  CMake compilation instructions for DiffElec tests                             #
# ------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------- #
# --               Configuration and Dependencies                -- #
# ----------------------------------------------------------------- #
# CMake version
cmake_minimum_required(VERSION 3.1)

# -- Dependency

find_package(Boost 1.48.0 REQUIRED)
include_directories(${BOOST_INCLUDE_DIRS})
set(LIBS ${LIBS} ${BOOST_LIBRARIES})

# -- Executables for tests
add_executable("main" "main.cpp")
target_link_libraries("main" DiffElec ${LIBS})
