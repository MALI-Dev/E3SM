
.SUFFIXES: .F .o .cpp
.PHONY: mode_forward shared analysis_members

FW = ../../../mpas-framework/src
SHARED_INCLUDES  = -I$(FW)/framework -I$(FW)/external/esmf_time_f90 -I$(FW)/operators
SHARED_INCLUDES += -I$(PWD)/shared -I$(PWD)/analysis_members -I$(PWD)/mode_forward
ifeq "$(SLM)" "true"
SHARED_INCLUDES += -I$(PWD)/SeaLevelModel
endif
ifeq "$(FASTISOSTASY)" "true"
SHARED_INCLUDES += -I$(PWD)/FastIsostasy/libisostasy/include
SHARED_INCLUDES += -I$(PWD)/fesm-utils/utils/include-serial
endif

FASTISOSTASY_FESM_MACHINE ?= linux
FASTISOSTASY_FESM_COMPILER ?= gfortran
FASTISOSTASY_CONFIG ?= config/foehn_gfortran
FASTISOSTASY_NC_INCLUDES ?= $(FCINCLUDES)
FASTISOSTASY_NC_LIBS ?= $(LIBS)

all: core_landice

shared: 
	(cd shared; $(MAKE) FCINCLUDES="$(FCINCLUDES) $(SHARED_INCLUDES)")

analysis_members: shared
	(cd analysis_members; $(MAKE) FCINCLUDES="$(FCINCLUDES) $(SHARED_INCLUDES)")


mode_forward: shared analysis_members sea_level_model fastisostasy_lib
	(cd mode_forward; $(MAKE) FCINCLUDES="$(FCINCLUDES) $(SHARED_INCLUDES)")

sea_level_model:
ifeq "$(SLM)" "true"
	if [ -e SeaLevelModel/.git ]; then \
		(cd SeaLevelModel; $(MAKE) all FC="$(FC)" FCFLAGS="$(FFLAGS)" FINCLUDES="$(FINCLUDES)") \
	else \
                (echo "Missing mpas-ocean/src/SeaLevelModel/.git, did you forget to 'git submodule update --init --recursive' ?"; exit 1) \
	fi
endif

fastisostasy_lib:
ifeq "$(FASTISOSTASY)" "true"
	if [ -e FastIsostasy/.git ]; then \
		if [ -e fesm-utils/.git ]; then \
			(cd fesm-utils; python3 build.py -m $(FASTISOSTASY_FESM_MACHINE) -c $(FASTISOSTASY_FESM_COMPILER) --component fftw --variant serial); \
			(cd fesm-utils; python3 build.py -m $(FASTISOSTASY_FESM_MACHINE) -c $(FASTISOSTASY_FESM_COMPILER) --component utils --variant serial); \
			(cd FastIsostasy; python3 config.py $(FASTISOSTASY_CONFIG)); \
			(cd FastIsostasy; $(MAKE) isostasy-static \
				FC="$(FC)" \
				INC_NC="$(FASTISOSTASY_NC_INCLUDES)" \
				LIB_NC="$(FASTISOSTASY_NC_LIBS)" \
				FESMUTILSROOT=../fesm-utils/utils \
				FFTWROOT=../fesm-utils/fftw-serial); \
		else \
			(echo "Missing components/mpas-albany-landice/src/fesm-utils/.git, did you forget to add or init the submodule?"; exit 1); \
		fi \
	else \
		(echo "Missing components/mpas-albany-landice/src/FastIsostasy/.git, did you forget to 'git submodule update --init --recursive'?"; exit 1); \
	fi
endif

core_landice: mode_forward shared analysis_members sea_level_model fastisostasy_lib
	ar -ru libdycore.a `find . -type f -name "*.o"`

core_input_gen:
	if [ ! -e default_inputs ]; then  mkdir default_inputs; fi
	(cd default_inputs; $(NL_GEN) ../Registry_processed.xml namelist.landice )
	(cd default_inputs; $(ST_GEN) ../Registry_processed.xml streams.landice stream_list.landice. listed )

core_reg:
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) Registry.xml > Registry_processed.xml

gen_includes:
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) Registry.xml > Registry_processed.xml
	(if [ ! -d inc ]; then mkdir -p inc; fi) # To generate *.inc files
	(cd inc; $(REG_PARSE) < ../Registry_processed.xml )

post_build:
	if [ ! -e $(ROOT_DIR)/default_inputs ]; then mkdir $(ROOT_DIR)/default_inputs; fi
	cp default_inputs/* $(ROOT_DIR)/default_inputs/.
	( cd $(ROOT_DIR)/default_inputs; for FILE in `ls -1`; do if [ ! -e ../$$FILE ]; then cp $$FILE ../.; fi; done )

clean:
	$(RM) *.o *.mod *.f90 libdycore.a
	$(RM) Registry_processed.xml
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i
	$(RM) -r default_inputs
	(cd shared; $(MAKE) clean)
	(cd mode_forward; $(MAKE) clean)
	(cd analysis_members; $(MAKE) clean)
	if [ -e SeaLevelModel/.git ]; then \
		(cd SeaLevelModel; $(MAKE) clean) \
	fi
	if [ -e FastIsostasy/.git ] && [ -e FastIsostasy/Makefile ]; then \
		(cd FastIsostasy; $(MAKE) clean) \
	fi
	if [ -e fesm-utils/.git ] && [ -e fesm-utils/utils/Makefile ]; then \
		(cd fesm-utils/utils; $(MAKE) clean) \
	fi
