set(CFLAGS "-mcmodel=medium")
if (compile_threaded)
  string(APPEND CFLAGS " -fopenmp")
endif()
if (DEBUG)
  string(APPEND CFLAGS " -g -Wall -fbacktrace -fcheck=bounds -ffpe-trap=invalid,zero,overflow")
endif()
if (NOT DEBUG)
  string(APPEND CFLAGS " -O")
endif()
if (COMP_NAME STREQUAL csm_share)
  string(APPEND CFLAGS " -std=c99")
endif()
set(CXXFLAGS "-std=c++14")
if (compile_threaded)
  string(APPEND CXXFLAGS " -fopenmp")
endif()
if (DEBUG)
  string(APPEND CXXFLAGS " -g -Wall -fbacktrace")
endif()
if (NOT DEBUG)
  string(APPEND CXXFLAGS " -O")
endif()
string(APPEND CPPDEFS " -DFORTRANUNDERSCORE -DNO_R16 -DCPRGNU")
if (DEBUG)
  string(APPEND CPPDEFS " -DYAKL_DEBUG")
endif()
set(CXX_LINKER "FORTRAN")
set(FC_AUTO_R8 "-fdefault-real-8")
set(FFLAGS "-mcmodel=medium -fconvert=big-endian -ffree-line-length-none -ffixed-line-length-none")
if (compile_threaded)
  string(APPEND FFLAGS " -fopenmp")
endif()
if (DEBUG)
  string(APPEND FFLAGS " -g -Wall -fbacktrace -fcheck=bounds -ffpe-trap=zero,overflow")
endif()
if (NOT DEBUG)
  string(APPEND FFLAGS " -O")
endif()
set(FFLAGS_NOOPT "-O0")
set(FIXEDFLAGS "-ffixed-form")
set(FREEFLAGS "-ffree-form")
set(HAS_F2008_CONTIGUOUS "FALSE")
if (compile_threaded)
  string(APPEND LDFLAGS " -fopenmp")
endif()
set(MPICC "mpicc")
set(MPICXX "mpicxx")
set(MPIFC "mpif90")
set(SCC "gcc")
set(SCXX "g++")
set(SFC "gfortran")
set(SUPPORTS_CXX "TRUE")
