# Load all kokkos settings from Ekat's mach file
set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)
include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)

# Enable Broadwell arch in Kokkos
option(Kokkos_ARCH_BDW "" ON)

set(SCREAM_MPIRUN_EXE "srun" CACHE STRING "")
set(SCREAM_MPI_NP_FLAG "-n" CACHE STRING "")
set(SCREAM_MPI_EXTRA_ARGS "" CACHE STRING "")
set(SCREAM_INPUT_ROOT "/usr/gdata/climdat/ccsm3data/inputdata" CACHE STRING "")

# Mach-specific flags
set(SCREAM_CXX_FLAGS "-w -cxxlib=/usr/tce/packages/gcc/gcc-8.3.1/rh" CACHE STRING "")
set(SCREAM_LD_FLAGS "-L/usr/tce/packages/gcc/gcc-8.3.1/rh/lib/gcc/x86_64-redhat-linux/8/ -mkl" CACHE STRING "")
