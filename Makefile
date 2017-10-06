# Generic Makefile for GEMMA
#
# Supported platforms
#
#       Unix / Linux               	LNX (default)
#       Mac                        	MAC
#
# Compilation options
#       static compilation    		FORCE_STATIC
#
# Examples:
#
#    Make GEMMA on Linux with OPENBLAS support:
#
#      make WITH_OPENBLAS=1
#
#    Disable debug info and checks (slightly faster release mode)
#
#      make WITH_OPENBLAS=1 DEBUG=
#
#    Force static compilation
#
#      make FORCE_STATIC=1
#
#    Run tests with
#
#      make check
#
#    Run quick (development) tests with
#
#      make fast-check
#
#    Run full (lengthy) tests with
#
#      make check-all
#
#    See also the INSTALL.md document in the source tree at
#
#      https://github.com/genetics-statistics/GEMMA/blob/master/INSTALL.md

# Set this variable to either LNX or MAC
SYS                    = LNX # LNX|MAC (Linux is the default)
# Leave blank after "=" to disable; put "= 1" to enable
DIST_NAME              = gemma-0.97.2
DEBUG                  = 1   # DEBUG mode, set DEBUG=0 for a release
SHOW_COMPILER_WARNINGS =
WITH_LAPACK            = 1
WITH_OPENBLAS          =     # Defaults to LAPACK - OPENBLAS may be faster
FORCE_STATIC           =     # Static linking of libraries
GCC_FLAGS              = -O3 # extra flags -Wl,--allow-multiple-definition
TRAVIS_CI              =     # used by TRAVIS for testing
EIGEN_INCLUDE_PATH=/usr/include/eigen3

# --------------------------------------------------------------------
# Edit below this line with caution
# --------------------------------------------------------------------

BIN_DIR  = ./bin
SRC_DIR  = ./src
TEST_SRC_DIR  = ./test/src

ifdef CXX # CXX defined in environment
  CPP = $(CXX)
  CC = $(CXX)
else
  CPP = g++
endif

ifdef OPENBLAS
  WITH_LAPACK =  # OPENBLAS usually includes LAPACK
endif

ifdef DEBUG
  CPPFLAGS = -g $(GCC_FLAGS) -std=gnu++11 -isystem/$(EIGEN_INCLUDE_PATH) -Icontrib/catch-1.9.7 -Isrc
else
  # release mode
  CPPFLAGS = -DNDEBUG $(GCC_FLAGS) -std=gnu++11 -isystem/$(EIGEN_INCLUDE_PATH) -Icontrib/catch-1.9.7 -Isrc
endif

ifdef SHOW_COMPILER_WARNINGS
  CPPFLAGS += -Wall
endif

ifndef FORCE_STATIC
  LIBS = -lgsl -lgslcblas -pthread -lz
else
  ifndef TRAVIS_CI # Travis static compile we cheat a little
    CPPFLAGS += -static
  endif
endif

OUTPUT = $(BIN_DIR)/gemma

SOURCES = $(SRC_DIR)/main.cpp

HDR =

# Detailed libary paths, D for dynamic and S for static

LIBS_LNX_D_LAPACK = -llapack
LIBS_LNX_D_BLAS = -lblas
LIBS_LNX_D_OPENBLAS = -lopenblas
LIBS_MAC_D_LAPACK = -framework Veclib
# LIBS_LNX_S_LAPACK = /usr/lib/libgsl.a  /usr/lib/libgslcblas.a /usr/lib/lapack/liblapack.a -lz
LIBS_LNX_S_LAPACK = /usr/lib/lapack/liblapack.a -lgfortran  /usr/lib/atlas-base/libatlas.a /usr/lib/libblas/libblas.a -Wl,--allow-multiple-definition

ifdef WITH_LAPACK
  ifeq ($(SYS), MAC)
    LIBS += $(LIBS_MAC_D_LAPACK)
  else
    ifndef FORCE_STATIC
      ifdef WITH_OPENBLAS
        LIBS += $(LIBS_LNX_D_OPENBLAS)
      else
        LIBS += $(LIBS_LNX_D_BLAS)
      endif
      LIBS += $(LIBS_LNX_D_LAPACK)
    else
      LIBS += $(LIBS_LNX_S_LAPACK)
    endif
  endif
endif

HDR          = $(wildcard src/*.h)
SOURCES      = $(wildcard src/*.cpp)

# all
OBJS = $(SOURCES:.cpp=.o)

all: $(OUTPUT)

$(OUTPUT): $(OBJS)
	$(CPP) $(CPPFLAGS) $(OBJS) $(LIBS) -o $(OUTPUT)

$(OBJS)  : $(HDR)

.cpp.o:
	$(CPP) $(CPPFLAGS) $(HEADERS) -c $*.cpp -o $*.o
.SUFFIXES : .cpp .c .o $(SUFFIXES)

unittests: all contrib/catch-1.9.7/catch.hpp $(TEST_SRC_DIR)/unittests-main.o $(TEST_SRC_DIR)/unittests-math.o
	$(CPP) $(CPPFLAGS) $(TEST_SRC_DIR)/unittests-main.o  $(TEST_SRC_DIR)/unittests-math.o $(filter-out src/main.o, $(OBJS)) $(LIBS) -o ./bin/unittests-gemma
	./bin/unittests-gemma

fast-check: all unittests
	rm -vf test/output/*
	cd test && ./dev_test_suite.sh | tee ../dev_test.log
	grep -q 'success rate: 100%' dev_test.log

slow-check: all
	rm -vf test/output/*
	cd test && ./test_suite.sh | tee ../test.log
	grep -q 'success rate: 100%' test.log

lengthy-check: all
	rm -vf test/output/*
	cd test && ./lengthy_test_suite.sh | tee ../lengthy_test.log
	grep -q 'success rate: 100%' lengthy_test.log

check: fast-check slow-check

check-all: check lengthy-check

clean:
	rm -vf $(SRC_DIR)/*.o
	rm -vf $(SRC_DIR)/*~
	rm -vf $(TEST_SRC_DIR)/*.o
	rm -vf $(OUTPUT)
	rm -vf ./bin/unittests-gemma

DIST_COMMON = COPYING.txt README.txt Makefile
DIST_SUBDIRS = src doc example bin

tar:
	mkdir -p ./$(DIST_NAME)
	cp $(DIST_COMMON) ./$(DIST_NAME)/
	cp -r $(DIST_SUBDIRS) ./$(DIST_NAME)/
	tar cvzf $(DIST_NAME).tar.gz ./$(DIST_NAME)/
	rm -r ./$(DIST_NAME)
