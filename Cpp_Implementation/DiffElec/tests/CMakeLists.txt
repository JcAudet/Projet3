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

# -- Dependency (Google Test)
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GTEST_LIBRARIES})


# -- Executables for tests
add_executable("UnitTest-Domain" "UnitTest-Domain.cpp")
target_link_libraries("UnitTest-Domain" DiffElec ${LIBS})
target_link_libraries("UnitTest-Domain" ${GTEST_LIBRARIES} pthread)
add_test(NAME "UnitTest-Domain" COMMAND $<TARGET_FILE:UnitTest-Domain>)
