set(CFLAGS "-O3 -qstrict")
string(APPEND CFLAGS " -qtune=440 -qarch=440d")
set(CONFIG_ARGS "--build=powerpc-bgp-linux --host=powerpc64-suse-linux")
string(APPEND CPPDEFS " -DLINUX -DnoI8")
string(APPEND FFLAGS " -qtune=440 -qarch=440d")
if (NOT DEBUG)
  string(APPEND FFLAGS " -O3 -qstrict -Q")
endif()
if (DEBUG)
  string(APPEND FFLAGS " -qinitauto=FF911299 -qflttrap=ov:zero:inv:en")
endif()
string(APPEND FFLAGS " -qextname=flush")
set(LDFLAGS "-Wl,--relax -Wl,--allow-multiple-definition")
set(MLIBS "-L/bgl/BlueLight/ppcfloor/bglsys/lib -lmpich.rts -lmsglayer.rts -lrts.rts -ldevices.rts")
set(MPICC "blrts_xlc")
set(MPIFC "blrts_xlf2003")
set(MPI_LIB_NAME "mpich.rts")
set(MPI_PATH "/bgl/BlueLight/ppcfloor/bglsys")
set(SCC "blrts_xlc")
set(SFC "blrts_xlf2003")
