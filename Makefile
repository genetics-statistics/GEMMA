#Makefile

# Supported platforms
#       Unix / Linux               	LNX
#       Mac                        	MAC
# Compilation options
#       32-bit binary        		FORCE_32BIT
#       dynamic compilation    		FORCE_DYNAMIC

# Set this variable to either LNX or MAC
SYS = LNX
# Leave blank after "=" to disable; put "= 1" to enable
SHOW_COMPILER_WARNINGS   =
WITH_LAPACK     = 1
WITH_OPENBLAS   =
NO_INTEL_COMPAT =
FORCE_32BIT     =
FORCE_DYNAMIC   =
DIST_NAME       = gemma-0.97

# --------------------------------------------------------------------
# Edit below this line with caution
# --------------------------------------------------------------------

EIGEN_INCLUDE_PATH=/usr/include/eigen3

BIN_DIR  = ./bin

SRC_DIR  = ./src

CPP = g++

ifdef DEBUG
  # development mode
  CPPFLAGS = -g -O3 -std=gnu++11 -isystem/$(EIGEN_INCLUDE_PATH)
else
  # release
  CPPFLAGS = -DNDEBUG -O3 -std=gnu++11 -isystem/$(EIGEN_INCLUDE_PATH)
endif

ifdef SHOW_COMPILER_WARNINGS
  CPPFLAGS += -Wall
endif

ifdef FORCE_DYNAMIC
  LIBS = -lgsl -lgslcblas -pthread -lz
else
  LIBS = -lgsl -lgslcblas -pthread -lz
endif

OUTPUT = $(BIN_DIR)/gemma

SOURCES = $(SRC_DIR)/main.cpp

HDR =

# Detailed libary paths, D for dynamic and S for static

LIBS_LNX_D_LAPACK = -llapack
LIBS_LNX_D_BLAS = -lblas
LIBS_LNX_D_OPENBLAS = -lopenblas
LIBS_MAC_D_LAPACK = -framework Veclib
LIBS_LNX_S_LAPACK = /usr/lib/lapack/liblapack.a -lgfortran  /usr/lib/atlas-base/libatlas.a /usr/lib/libblas/libblas.a -Wl,--allow-multiple-definition

SOURCES += $(SRC_DIR)/param.cpp $(SRC_DIR)/gemma.cpp $(SRC_DIR)/io.cpp $(SRC_DIR)/lm.cpp $(SRC_DIR)/lmm.cpp $(SRC_DIR)/vc.cpp $(SRC_DIR)/mvlmm.cpp $(SRC_DIR)/bslmm.cpp $(SRC_DIR)/prdt.cpp $(SRC_DIR)/mathfunc.cpp $(SRC_DIR)/gzstream.cpp $(SRC_DIR)/eigenlib.cpp $(SRC_DIR)/ldr.cpp $(SRC_DIR)/bslmmdap.cpp $(SRC_DIR)/logistic.cpp $(SRC_DIR)/varcov.cpp
HDR += $(SRC_DIR)/param.h $(SRC_DIR)/gemma.h $(SRC_DIR)/io.h $(SRC_DIR)/lm.h $(SRC_DIR)/lmm.h $(SRC_DIR)/vc.h $(SRC_DIR)/mvlmm.h $(SRC_DIR)/bslmm.h $(SRC_DIR)/prdt.h $(SRC_DIR)/mathfunc.h $(SRC_DIR)/gzstream.h $(SRC_DIR)/eigenlib.h

ifdef WITH_LAPACK
  OBJS += $(SRC_DIR)/lapack.o
ifeq ($(SYS), MAC)
  LIBS += $(LIBS_MAC_D_LAPACK)
else
  ifdef FORCE_DYNAMIC
    ifdef WITH_OPENBLAS
      LIBS += $(LIBS_LNX_D_OPENBLAS)
    else
      LIBS += $(LIBS_LNX_D_BLAS)
    endif
    LIBS += $(LIBS_LNX_D_LAPACK)
  endif
endif
  SOURCES += $(SRC_DIR)/lapack.cpp
  HDR += $(SRC_DIR)/lapack.h
endif

ifdef NO_INTEL_COMPAT
  else
  ifdef FORCE_32BIT
    CPPFLAGS += -m32
  else
    CPPFLAGS += -m64
  endif
endif

ifdef FORCE_DYNAMIC
else
  CPPFLAGS += -static
endif


# all
OBJS = $(SOURCES:.cpp=.o)

all: $(OUTPUT)

$(OUTPUT): $(OBJS)
	$(CPP) $(CPPFLAGS) $(OBJS) $(LIBS) -o $(OUTPUT)

$(OBJS) : $(HDR)

.cpp.o:
	$(CPP) $(CPPFLAGS) $(HEADERS) -c $*.cpp -o $*.o
.SUFFIXES : .cpp .c .o $(SUFFIXES)

fast-check: all
	cd test && ./dev_test_suite.sh | tee ../dev_test.log
	grep -q 'success rate: 100%' dev_test.log

slow-check: all
	cd test && ./test_suite.sh | tee ../test.log
	grep -q 'success rate: 100%' test.log

check: fast-check slow-check

clean:
	rm -rf ${SRC_DIR}/*.o ${SRC_DIR}/*~ *~ $(OUTPUT)
	rm test/output/*

DIST_COMMON = COPYING.txt README.txt Makefile
DIST_SUBDIRS = src doc example bin

tar:
	mkdir -p ./$(DIST_NAME)
	cp $(DIST_COMMON) ./$(DIST_NAME)/
	cp -r $(DIST_SUBDIRS) ./$(DIST_NAME)/
	tar cvzf $(DIST_NAME).tar.gz ./$(DIST_NAME)/
	rm -r ./$(DIST_NAME)
