KOKKOS_PATH = ../../../kokkos
KOKKOS_DEVICES = "Cuda"
EXE_NAME = "02_Exercise"

SRC = $(wildcard *.cpp)
BUILDDIR = ../../../../build/simtbx/kokkos
LIBDIR = ../../../../build/lib

LIB_TARGET = libsimtbx_kokkos.so
TARGET = simtbx_kokkos_ext.so

CCTBXDIR = $(shell pwd)/../../../../

default: build
	echo "Start Build"


ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
KOKKOS_ARCH = "Volta70"
KOKKOS_CUDA_OPTIONS = "enable_lambda,force_uvm"
else
CXX = g++
KOKKOS_ARCH = "HSW"
endif

CXXFLAGS = -O3 -fPIC -DCUDAREAL=double
LINK = g++
LINKFLAGS = -s -std=c++11 -fopenmp
LDFLAGS = -Llib -L$(CUDA_HOME)/lib64 
LDFLAGS += -L$(LIBDIR)
LDFLAGS += -L$(CCTBXDIR)/modules/lib
LDFLAGS += -L$(CCTBXDIR)/modules/cctbx_project/lib
LDFLAGS += -L$(CCTBXDIR)/modules/cctbx_project/simtbx/kokkos
LDLIBS = -lkokkos -ldl -lcudart -lcuda -lboost_thread -lboost_system -lboost_python -lm

EXTRA_LDLIBS = -lsimtbx_kokkos -lomptbx -lscitbx_boost_python -lcctbx 

INCLUDES = -I../../../../modules/boost
INCLUDES += -I../../../../build/include
INCLUDES += -I../../../../modules/cctbx_project
INCLUDES += -I../../../../modules
INCLUDES += -I../../../../opt/mamba/envs/psana_env/include/python3.7m

DEPFLAGS = -M

OBJ = $(SRC:%.cpp=$(BUILDDIR)/%.o)

include $(KOKKOS_PATH)/Makefile.kokkos

build: dir $(LIBDIR)/$(TARGET)

dir: 
	mkdir -p $(BUILDDIR)

# we really only need 'kokkos_ext.so', not whole OBJ. Easier way?
$(LIBDIR)/$(TARGET): $(LIBDIR)/$(LIB_TARGET) $(KOKKOS_LINK_DEPENDS)
	$(LINK) -o $@ -shared $(LINKFLAGS) $(OBJ) $(LDFLAGS) $(EXTRA_LDLIBS) $(LDLIBS) 

clean: kokkos-clean 
	rm -f $(BUILDDIR)/*.o $(BUILDDIR)/*.cuda $(BUILDDIR)/*.host

# Compilation rules

$(LIBDIR)/$(LIB_TARGET): $(OBJ) $(KOKKOS_LINK_DEPENDS)
	$(CXX) -o $@ -shared $(OBJ) $(LDFLAGS) $(LDLIBS) 

$(OBJ): $(BUILDDIR)/%.o: %.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(INCLUDES) $(CXXFLAGS) $(EXTRA_INC) -c $< -o $@
