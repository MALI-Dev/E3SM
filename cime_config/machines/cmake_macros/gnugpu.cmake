string(APPEND CFLAGS " -mcmodel=medium")
string(APPEND FFLAGS " -mcmodel=medium -fconvert=big-endian -ffree-line-length-none -ffixed-line-length-none")
if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
   string(APPEND FFLAGS " -fallow-argument-mismatch")
endif()
string(APPEND CPPDEFS " -DFORTRANUNDERSCORE -DNO_R16 -DCPRGNU")
if (compile_threaded)
  string(APPEND CFLAGS " -fopenmp")
  string(APPEND CXXFLAGS " -fopenmp")
  string(APPEND FFLAGS " -fopenmp")
  string(APPEND LDFLAGS " -fopenmp")
endif()
if (DEBUG)
  string(APPEND CFLAGS   " -O0 -g -Wall -fbacktrace -fcheck=bounds -ffpe-trap=invalid,zero,overflow")
  string(APPEND FFLAGS   " -O0 -g -Wall -fbacktrace -fcheck=bounds -ffpe-trap=zero,overflow")
  string(APPEND CXXFLAGS " -O0 -g -Wall -fbacktrace")
  string(APPEND CPPDEFS  " -DYAKL_DEBUG")
endif()
if (NOT DEBUG)
  string(APPEND CFLAGS " -O")
  string(APPEND CXXFLAGS " -O")
  string(APPEND FFLAGS " -O")
endif()
if (COMP_NAME STREQUAL csm_share)
  string(APPEND CFLAGS " -std=c99")
endif()
set(SUPPORTS_CXX "TRUE")
string(APPEND CXX_LIBS " -lstdc++")
string(APPEND FFLAGS_NOOPT " -O0")
string(APPEND FIXEDFLAGS " -ffixed-form")
string(APPEND FREEFLAGS " -ffree-form")
set(HAS_F2008_CONTIGUOUS "FALSE")
set(MPICC "mpicc")
set(MPICXX "mpicxx")
set(MPIFC "mpif90")
set(SCC "gcc")
set(SCXX "g++")
set(SFC "gfortran")
