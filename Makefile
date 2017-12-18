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
#    Make GEMMA on Linux without OPENBLAS support:
#
#      make WITH_OPENBLAS=
#
#    Disable debug info and checks (slightly faster release mode)
#
#      make DEBUG=
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

GEMMA_VERSION = $(shell cat ./VERSION)

# Set this variable to either LNX or MAC
SYS                    = LNX # LNX|MAC (Linux is the default)
# Leave blank after "=" to disable; put "= 1" to enable
DIST_NAME              = gemma-$(GEMMA_VERSION)
DEBUG                  = 1                # DEBUG mode, set DEBUG=0 for a release
SHOW_COMPILER_WARNINGS =
WITH_LAPACK            = 1
WITH_OPENBLAS          = 1                # Without OpenBlas uses LAPACK
OPENBLAS_LEGACY        =                  # Using older OpenBlas
FORCE_STATIC           =                  # Static linking of libraries
GCC_FLAGS              = -Wall -O3 -std=gnu++11 # extra flags -Wl,--allow-multiple-definition
TRAVIS_CI              =                  # used by TRAVIS for testing
EIGEN_INCLUDE_PATH     = /usr/include/eigen3
OPENBLAS_INCLUDE_PATH  = /usr/local/opt/openblas/include

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

ifeq ($(CPP), clang++)
  # macOS Homebrew settings (as used on Travis-CI)
  GCC_FLAGS=-O3 -std=c++11 -stdlib=libc++ -isystem/$(OPENBLAS_INCLUDE_PATH) -isystem//usr/local/include/eigen3 -Wl,-L/usr/local/opt/openblas/lib
endif

ifdef WITH_OPENBLAS
  OPENBLAS=1
  # WITH_LAPACK =  # OPENBLAS usually includes LAPACK
  CPPFLAGS += -DOPENBLAS -isystem/$(OPENBLAS_INCLUDE_PATH)
  ifdef OPENBLAS_LEGACY
    # Legacy version (mostly for Travis-CI)
    CPPFLAGS += -DOPENBLAS_LEGACY
  endif
endif

ifdef DEBUG
  CPPFLAGS += -g $(GCC_FLAGS) -isystem/$(EIGEN_INCLUDE_PATH) -Icontrib/catch-1.9.7 -Isrc
else
  # release mode
  CPPFLAGS += -DNDEBUG $(GCC_FLAGS) -isystem/$(EIGEN_INCLUDE_PATH) -Icontrib/catch-1.9.7 -Isrc
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

.PHONY: all

OUTPUT = $(BIN_DIR)/gemma

# Detailed libary paths, D for dynamic and S for static

# LIBS_LNX_D_LAPACK = -llapack
# LIBS_LNX_D_BLAS = -lblas
LIBS_LNX_D_OPENBLAS = -lopenblas
LIBS_MAC_D_LAPACK = -framework Accelerate
# LIBS_LNX_S_LAPACK = /usr/lib/libgsl.a  /usr/lib/libgslcblas.a /usr/lib/lapack/liblapack.a -lz
# LIBS_LNX_S_LAPACK = /usr/lib/lapack/liblapack.a -lgfortran  /usr/lib/atlas-base/libatlas.a /usr/lib/libblas/libblas.a -Wl,--allow-multiple-definition

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

HDR          = $(wildcard src/*.h) ./src/version.h
SOURCES      = $(wildcard src/*.cpp)

# all
OBJS = $(SOURCES:.cpp=.o)

all: $(OUTPUT)

./src/version.h:
	./scripts/gen_version_info.sh > src/version.h

$(OUTPUT): $(OBJS)
	$(CPP) $(CPPFLAGS) $(OBJS) $(LIBS) -o $(OUTPUT)

$(OBJS): $(HDR)

# .cpp.o:
# 	$(CPP) $(CPPFLAGS) -c $*.cpp -o $*.o

.SUFFIXES : .cpp .c .o $(SUFFIXES)

./bin/unittests-gemma: contrib/catch-1.9.7/catch.hpp $(TEST_SRC_DIR)/unittests-main.o $(TEST_SRC_DIR)/unittests-math.o $(OBJS)
	$(CPP) $(CPPFLAGS) $(TEST_SRC_DIR)/unittests-main.o  $(TEST_SRC_DIR)/unittests-math.o $(filter-out src/main.o, $(OBJS)) $(LIBS) -o ./bin/unittests-gemma

unittests: ./bin/unittests-gemma
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
	rm $(SRC_DIR)/version.h
	rm -vf $(SRC_DIR)/*.o
	rm -vf $(SRC_DIR)/*~
	rm -vf $(TEST_SRC_DIR)/*.o
	rm -vf $(OUTPUT)
	rm -vf ./bin/unittests-gemma

DIST_COMMON = *.md LICENSE VERSION Makefile
DIST_SUBDIRS = src doc example bin

tar: version all
	@echo "Creating $(DIST_NAME)"
	mkdir -p ./$(DIST_NAME)
	cp $(DIST_COMMON) ./$(DIST_NAME)/
	cp -r $(DIST_SUBDIRS) ./$(DIST_NAME)/
	tar cvzf $(DIST_NAME).tar.gz ./$(DIST_NAME)/
	rm -r ./$(DIST_NAME)
