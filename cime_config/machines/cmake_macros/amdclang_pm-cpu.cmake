if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_SLASHPROC -DHAVE_GETTIMEOFDAY")
endif()
if (NOT DEBUG)
  string(APPEND CFLAGS " -O2 -g")
  string(APPEND CXXFLAGS " -O2 -g")
  string(APPEND FFLAGS " -O2")
endif()
#string(APPEND FFLAGS " -march=znver3")
set(SCC "clang")
set(SCXX "clang++")
set(SFC "flang")

if (COMP_NAME STREQUAL cism)
  string(APPEND CMAKE_OPTS " -D CISM_GNU=ON")
endif()

string(APPEND FFLAGS " -Mflushz ")
string(APPEND FIXEDFLAGS " -Mfixed")
string(APPEND FREEFLAGS " -Mfreeform")
if (compile_threaded)
  string(APPEND FFLAGS " -mp")
  string(APPEND LDFLAGS " -mp")
endif()
