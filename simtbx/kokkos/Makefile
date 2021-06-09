KOKKOS_PATH = ../../../kokkos
KOKKOS_DEVICES = "Cuda"
EXE_NAME = "02_Exercise"

SRC = $(wildcard *.cpp)
BUILDDIR = ../../../../build/simtbx/kokkos

default: build
	echo "Start Build"

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
EXE = ${EXE_NAME}.cuda
KOKKOS_ARCH = "Volta70"
KOKKOS_CUDA_OPTIONS = "enable_lambda,force_uvm"
else
CXX = g++
EXE = ${EXE_NAME}.host
KOKKOS_ARCH = "BDW"
endif

CXXFLAGS = -O3 -fPIC -DCUDAREAL=double
LINK = ${CXX}
LINKFLAGS =

INCLUDES = -I../../../../modules/boost
INCLUDES += -I../../../../build/include
INCLUDES += -I../../../../modules/cctbx_project
INCLUDES += -I../../../../modules
INCLUDES += -I../../../../opt/mamba/envs/psana_env/include/python3.7m

DEPFLAGS = -M

OBJ = $(SRC:%.cpp=$(BUILDDIR)/%.o)
LIB =

include $(KOKKOS_PATH)/Makefile.kokkos

build: dir $(BUILDDIR)/$(EXE)

dir: 
	mkdir -p $(BUILDDIR)

$(BUILDDIR)/$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) $(LIB) -o $@

clean: kokkos-clean 
	rm -f $(BUILDDIR)/*.o $(BUILDDIR)/*.cuda $(BUILDDIR)/*.host

# Compilation rules

$(OBJ): $(BUILDDIR)/%.o: %.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(INCLUDES) $(CXXFLAGS) $(EXTRA_INC) -c $< -o $@
