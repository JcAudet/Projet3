# --------------------------------------------------------------------------- #
# Author:       Francois Fillion-Gourdeau             #
# Date:         2018-02-02                                                    #
# Description:  CMake compilation instructions.                               #
# ----------------------------------------------------------------------------#

# ----------------------------------------------------------------- #
# --                 Name and version of library                 -- #
# ----------------------------------------------------------------- #
project(KbPro)
set (KbPro_VERSION_MAJOR 0)
set (KbPro_VERSION_MINOR 1)
set (KbPro_VERSION_RELEASE 0)

# ----------------------------------------------------------------- #
# --               Configuration and Dependencies                -- #
# ----------------------------------------------------------------- #
# -- CMake version and installation directory.
cmake_minimum_required(VERSION 3.0.2)
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX /usr)
else()
    include_directories(${CMAKE_INSTALL_PREFIX}/include)
endif()
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

# -- Required dependencies: 

find_package(Boost 1.48.0 REQUIRED)
include_directories(${BOOST_INCLUDE_DIRS})
set(LIBS ${LIBS} ${BOOST_LIBRARIES})


#find_package(MPI REQUIRED)
#include_directories(${MPI_INCLUDE_PATH})
#link_libraries(${MPI_LIBRARIES})

find_package(OpenMP REQUIRED)
include_directories(${OpenMP_INCLUDE_PATH})
set(LIBS ${LIBS} ${OpenMP_LIBRARIES})
if(OPENMP_FOUND)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
#link_libraries(${OpenMP_LIBRARIES})

find_package(armadillo REQUIRED)
include_directories(${armadillo_INCLUDE_DIRS})
set(LIBS ${LIBS} ${armadillo_LIBRARIES})

find_package(FFTW3 REQUIRED)
include_directories(${FFTW3_INCLUDE_DIRS})
set (LIBS ${LIBS} ${FFTW_LIBRARIES})

#include_directories(${FFTW_MPI_INCLUDE_DIR})
#link_libraries(${FFTW_MPI_LIBRARIES})# -lfftw3_mpi) 

include_directories(${FFTW_OMP_INCLUDE_DIR})
link_libraries(${FFTW_OMP_LIBRARIES})# -lfftw3_omp) 

#find_package(HDF5 REQUIRED)
#include_directories(${HDF5_INCLUDE_DIRS})
#set (LIBS ${LIBS} ${HDF5_LIBRARIES})



#enable_language (Fortran)

# -- Compiler config (flags and such)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -march=native -std=c++11")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall -march=native")



#find_package(HDF5 REQUIRED)
#include_directories(${HDF5_INCLUDE_DIRS})
#set (LIBS ${LIBS} ${HDF5_LIBRARIES})

# ----------------------------------------------------------------- #
# --                  Compilation Instructions                   -- #
# ----------------------------------------------------------------- #
# -- Included files
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/include INC_LIST)

# -- Source files
aux_source_directory (${CMAKE_CURRENT_SOURCE_DIR}/src SRC_LIST)
aux_source_directory (${CMAKE_CURRENT_SOURCE_DIR}/src/fft SRC_LIST)
aux_source_directory (${CMAKE_CURRENT_SOURCE_DIR}/src/ComputeSelfEnergy SRC_LIST)
aux_source_directory (${CMAKE_CURRENT_SOURCE_DIR}/src/TwoTimesPropagation SRC_LIST)


# -- Build a shared library
# -- Build a shared library
add_library(${PROJECT_NAME} SHARED ${SRC_LIST})


# Define the properties of the library to
# be built (version and other stuff).
set_target_properties (
  ${PROJECT_NAME}
  PROPERTIES
    VERSION   ${KbPro_VERSION_MAJOR}.${KbPro_VERSION_MINOR}.${KbPro_VERSION_RELEASE}
    SOVERSION ${KbPro_VERSION_MAJOR}.${KbPro_VERSION_MINOR}.${KbPro_VERSION_RELEASE}
)

# Link dependencies
target_link_libraries(${PROJECT_NAME} ${LIBS})

# Install directories
install (TARGETS ${PROJECT_NAME} DESTINATION lib)
install (DIRECTORY include/ DESTINATION include)

# -- KbPro tests
include(CTest)
add_subdirectory(tests)

# -- KbPro drivers
#add_subdirectory(drivers/TimeIndependent/SingleComp)
#add_subdirectory(drivers/TimeDependent/SingleComp)
#add_subdirectory(drivers/TimeDependent/SingleCompRegCoulomb)

# -- Python library
#add_subdirectory(python_module)

# ----------------------------------------------------------------- #
# --                   Generate Documentation                    -- #
# ----------------------------------------------------------------- #
# add a target to generate API documentation with Doxygen
#find_package(Doxygen)
#if(DOXYGEN_FOUND)
#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
#add_custom_target(
#  doc ALL ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
#      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#      COMMENT "Generating API documentation with Doxygen" VERBATIM
#)
#endif(DOXYGEN_FOUND)
