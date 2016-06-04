#Makefile

# Supported platforms
#       Unix / Linux               	LNX
#       Mac                        	MAC
# Compilation options
#       link to LAPACK              WITH_LAPACK
#       32-bit binary        		FORCE_32BIT
#       dynamic compilation    		FORCE_DYNAMIC
#       float precision				FORCE_FLOAT

# Set this variable to either LNX or MAC
SYS = LNX
# Leave blank after "=" to disable; put "= 1" to enable
# Disable WITH_LAPACK option can slow computation speed significantly and is not recommended
# Disable WITH_ARPACK option only disable -apprx option in the software
WITH_LAPACK = 1
FORCE_32BIT = 
FORCE_DYNAMIC = 
FORCE_FLOAT = 
DIST_NAME = gemma-0.95alpha

# --------------------------------------------------------------------
# Edit below this line with caution
# --------------------------------------------------------------------


BIN_DIR  = ./bin

SRC_DIR  = ./src

CPP = g++

CPPFLAGS = -Wall -O3 -std=gnu++11

ifdef FORCE_DYNAMIC
LIBS = -lgsl -lgslcblas -lblas -pthread -lz
else
LIBS = -lgsl -lgslcblas -pthread -lz
endif

OUTPUT = $(BIN_DIR)/gemma

ifdef FORCE_FLOAT
OUTPUT = $(BIN_DIR)/gemmaf
endif

SOURCES = $(SRC_DIR)/main.cpp

HDR = 

# Detailed libary paths, D for dynamic and S for static

LIBS_LNX_D_LAPACK = -llapack
LIBS_MAC_D_LAPACK = -framework Accelerate
LIBS_LNX_S_LAPACK = /usr/lib/lapack/liblapack.a -lgfortran  /usr/lib/atlas-base/libatlas.a /usr/lib/libblas/libblas.a -Wl,--allow-multiple-definition 

# Options

ifdef FORCE_FLOAT
  CPPFLAGS += -DFORCE_FLOAT
  SOURCES += $(SRC_DIR)/param_float.cpp $(SRC_DIR)/gemma_float.cpp $(SRC_DIR)/io_float.cpp $(SRC_DIR)/lm_float.cpp $(SRC_DIR)/vc_float.cpp $(SRC_DIR)/lmm_float.cpp $(SRC_DIR)/mvlmm_float.cpp $(SRC_DIR)/bslmm_float.cpp $(SRC_DIR)/prdt_float.cpp $(SRC_DIR)/mathfunc_float.cpp $(SRC_DIR)/gzstream.cpp $(SRC_DIR)/eigenlib.cpp
  HDR += $(SRC_DIR)/param_float.h $(SRC_DIR)/gemma_float.h $(SRC_DIR)/io_float.h $(SRC_DIR)/lm_float.h $(SRC_DIR)/lmm_float.h $(SRC_DIR)/vc_float.h $(SRC_DIR)/mvlmm_float.h $(SRC_DIR)/bslmm_float.h $(SRC_DIR)/prdt_float.h $(SRC_DIR)/mathfunc_float.h $(SRC_DIR)/gzstream.h $(SRC_DIR)/eigenlib.h
else
  SOURCES += $(SRC_DIR)/param.cpp $(SRC_DIR)/gemma.cpp $(SRC_DIR)/io.cpp $(SRC_DIR)/lm.cpp $(SRC_DIR)/lmm.cpp $(SRC_DIR)/vc.cpp $(SRC_DIR)/mvlmm.cpp $(SRC_DIR)/bslmm.cpp $(SRC_DIR)/prdt.cpp $(SRC_DIR)/mathfunc.cpp $(SRC_DIR)/gzstream.cpp $(SRC_DIR)/eigenlib.cpp $(SRC_DIR)/utils.cpp
  HDR += $(SRC_DIR)/param.h $(SRC_DIR)/gemma.h $(SRC_DIR)/io.h $(SRC_DIR)/lm.h $(SRC_DIR)/lmm.h $(SRC_DIR)/vc.h $(SRC_DIR)/mvlmm.h $(SRC_DIR)/bslmm.h $(SRC_DIR)/prdt.h $(SRC_DIR)/mathfunc.h $(SRC_DIR)/gzstream.h $(SRC_DIR)/eigenlib.h $(SRC_DIR)/utils.h
endif

ifdef WITH_LAPACK
  OBJS += $(SRC_DIR)/lapack.o
  CPPFLAGS += -DWITH_LAPACK
ifeq ($(SYS), MAC)
  LIBS += $(LIBS_MAC_D_LAPACK)
else
ifdef FORCE_DYNAMIC
  LIBS += $(LIBS_LNX_D_LAPACK)
else
  LIBS += $(LIBS_LNX_S_LAPACK)
endif
endif
  SOURCES += $(SRC_DIR)/lapack.cpp
  HDR += $(SRC_DIR)/lapack.h
endif

ifdef FORCE_32BIT
  CPPFLAGS += -m32
else
  CPPFLAGS += -m64
endif

ifneq ($(SYS), MAC)
ifdef FORCE_DYNAMIC
else
  CPPFLAGS += -static
endif
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


clean:
	rm -rf ${SRC_DIR}/*.o ${SRC_DIR}/*~ *~ ${SRC_DIR}/*_float.* $(OUTPUT)

DIST_COMMON = COPYING.txt README.txt Makefile
DIST_SUBDIRS = src doc example bin

tar:
	mkdir -p ./$(DIST_NAME)
	cp $(DIST_COMMON) ./$(DIST_NAME)/
	cp -r $(DIST_SUBDIRS) ./$(DIST_NAME)/
	tar cvzf $(DIST_NAME).tar.gz ./$(DIST_NAME)/
	rm -r ./$(DIST_NAME)

