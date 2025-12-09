CXX       := clang++
MPICXX    := mpicxx
INCLUDES  := -Iutils/include -Isequential/include -Iparallel_openmp/include -Impi/include
LIBS      :=
LDFLAGS   :=

MPICXX_XL := mpixlC
CXXFLAGS_PWR8 := -O3 -qarch=pwr8 -qtune=pwr8 -qhot -qsimd=auto -std=c++11 $(INCLUDES)

MPICXX_XL_OMP := mpixlC
CXXFLAGS_HYBRID := -O3 -qarch=pwr8 -qtune=pwr8 -qhot -qsimd=auto -qsmp=omp -std=c++11 $(INCLUDES)

MPICXX_OPENMPI := mpicxx
CXXFLAGS_OPENMPI_HYBRID := -O3 -fopenmp -std=c++11 $(INCLUDES)

UTILS_SRC := utils/src/arg_parser.cpp utils/src/error_metrics.cpp
UTILS_OBJ := $(UTILS_SRC:.cpp=.o)

SEQUENTIAL_SRC := sequential/src/main.cpp sequential/src/solver.cpp
SEQUENTIAL_OBJ := $(filter-out sequential/src/main.o,$(SEQUENTIAL_SRC:.cpp=.o))

UTILS_OBJ_OMP := $(UTILS_SRC:.cpp=.omp.o)
SEQUENTIAL_OBJ_OMP := $(filter-out sequential/src/main.omp.o,$(SEQUENTIAL_SRC:.cpp=.omp.o))

MPI_SRC := mpi/src/main.cpp mpi/src/mpi_solver.cpp
MPI_OBJ := $(filter-out mpi/src/main.mpi.o,$(MPI_SRC:.cpp=.mpi.o))
UTILS_OBJ_MPI := $(UTILS_SRC:.cpp=.mpi.o)

MPI_OBJ_PWR8 := $(filter-out mpi/src/main.pwr8.o,$(MPI_SRC:.cpp=.pwr8.o))
UTILS_OBJ_PWR8 := $(UTILS_SRC:.cpp=.pwr8.o)

MPI_OBJ_HYBRID := $(filter-out mpi/src/main.hybrid.o,$(MPI_SRC:.cpp=.hybrid.o))
UTILS_OBJ_HYBRID := $(UTILS_SRC:.cpp=.hybrid.o)

MPI_OBJ_OPENMPI := $(filter-out mpi/src/main.openmpi.o,$(MPI_SRC:.cpp=.openmpi.o))
UTILS_OBJ_OPENMPI := $(UTILS_SRC:.cpp=.openmpi.o)

BUILD_DIR := build

CXXFLAGS_COMMON := -std=c++11 -Wall -Wextra $(INCLUDES)

CXXFLAGS_RELEASE := -O3

CXXFLAGS_PROFILE := -O3 -g -gline-tables-only -fno-omit-frame-pointer

CXXFLAGS_PROFILE_DEEP := $(CXXFLAGS_PROFILE) -fno-inline-functions

CXXFLAGS_DEBUG := -O0 -g3 -fno-omit-frame-pointer

BUILD ?= release
ifeq ($(BUILD),release)
  CXXFLAGS := $(CXXFLAGS_COMMON) $(CXXFLAGS_RELEASE)
else ifeq ($(BUILD),profile)
  CXXFLAGS := $(CXXFLAGS_COMMON) $(CXXFLAGS_PROFILE)
else ifeq ($(BUILD),profile-deep)
  CXXFLAGS := $(CXXFLAGS_COMMON) $(CXXFLAGS_PROFILE_DEEP)
else ifeq ($(BUILD),debug)
  CXXFLAGS := $(CXXFLAGS_COMMON) $(CXXFLAGS_DEBUG)
else
  $(error Unknown BUILD=$(BUILD))
endif

POSTLINK := @dsymutil $@ -o $@.dSYM >/dev/null 2>&1 || true

.PHONY: all clean sequential openmp mpi mpi_pwr8 mpi_hybrid_pwr8 mpi_hybrid_pwr8_openmpi \
        profile profile-deep debug openmp-profile openmp-profile-deep mpi-profile mpi-debug \
        trace trace-omp run run-omp run-mpi run-mpi-pwr8 run-mpi-hybrid run-mpi-hybrid-openmpi

all: sequential openmp mpi

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.omp.o: %.cpp
	$(CXX) $(CXXFLAGS) -fopenmp -c $< -o $@

%.mpi.o: %.cpp
	$(MPICXX) $(CXXFLAGS) -c $< -o $@

%.pwr8.o: %.cpp
	$(MPICXX_XL) $(CXXFLAGS_PWR8) -c $< -o $@

%.hybrid.o: %.cpp
	$(MPICXX_XL_OMP) $(CXXFLAGS_HYBRID) -c $< -o $@

%.openmpi.o: %.cpp
	$(MPICXX_OPENMPI) $(CXXFLAGS_OPENMPI_HYBRID) -c $< -o $@

sequential: $(BUILD_DIR)/sequential

$(BUILD_DIR)/sequential: sequential/src/main.cpp $(SEQUENTIAL_OBJ) $(UTILS_OBJ) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) sequential/src/main.cpp $(SEQUENTIAL_OBJ) $(UTILS_OBJ) \
		$(LDFLAGS) $(LIBS) -o $@
	$(POSTLINK)

profile:
	$(MAKE) sequential BUILD=profile
profile-deep:
	$(MAKE) sequential BUILD=profile-deep
debug:
	$(MAKE) sequential BUILD=debug

run: sequential
	./$(BUILD_DIR)/sequential $(ARGS)

trace: sequential
	xcrun xctrace record \
	  --template 'Time Profiler' \
	  --output $(BUILD_DIR)/timeprof.trace \
	  --launch -- ./$(BUILD_DIR)/sequential $(ARGS)

openmp: $(BUILD_DIR)/sequential_omp

$(BUILD_DIR)/sequential_omp: sequential/src/main.cpp $(SEQUENTIAL_OBJ_OMP) $(UTILS_OBJ_OMP) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -fopenmp sequential/src/main.cpp $(SEQUENTIAL_OBJ_OMP) $(UTILS_OBJ_OMP) \
		$(LDFLAGS) $(LIBS) -fopenmp -o $@
	$(POSTLINK)

openmp-profile:
	$(MAKE) openmp BUILD=profile
openmp-profile-deep:
	$(MAKE) openmp BUILD=profile-deep

run-omp: openmp
	./$(BUILD_DIR)/sequential_omp $(ARGS)

trace-omp: openmp
	xcrun xctrace record \
	  --template 'Time Profiler' \
	  --output $(BUILD_DIR)/timeprof_omp.trace \
	  --launch -- ./$(BUILD_DIR)/sequential_omp $(ARGS)

mpi: $(BUILD_DIR)/mpi_solver

$(BUILD_DIR)/mpi_solver: mpi/src/main.cpp $(MPI_OBJ) $(UTILS_OBJ_MPI) | $(BUILD_DIR)
	$(MPICXX) $(CXXFLAGS) mpi/src/main.cpp $(MPI_OBJ) $(UTILS_OBJ_MPI) \
		$(LDFLAGS) $(LIBS) -o $@
	$(POSTLINK)

mpi-profile:
	$(MAKE) mpi BUILD=profile
mpi-debug:
	$(MAKE) mpi BUILD=debug

run-mpi: mpi
	mpirun -np $(NP) ./$(BUILD_DIR)/mpi_solver $(ARGS)


mpi_pwr8: $(BUILD_DIR)/mpi_solver_pwr8

$(BUILD_DIR)/mpi_solver_pwr8: mpi/src/main.cpp $(MPI_OBJ_PWR8) $(UTILS_OBJ_PWR8) | $(BUILD_DIR)
	$(MPICXX_XL) $(CXXFLAGS_PWR8) mpi/src/main.cpp $(MPI_OBJ_PWR8) $(UTILS_OBJ_PWR8) \
		$(LDFLAGS) $(LIBS) -o $@

run-mpi-pwr8: mpi_pwr8
	mpirun -np $(NP) ./$(BUILD_DIR)/mpi_solver_pwr8 $(ARGS)

mpi_hybrid_pwr8: $(BUILD_DIR)/mpi_solver_hybrid

$(BUILD_DIR)/mpi_solver_hybrid: mpi/src/main.cpp $(MPI_OBJ_HYBRID) $(UTILS_OBJ_HYBRID) | $(BUILD_DIR)
	$(MPICXX_XL_OMP) $(CXXFLAGS_HYBRID) mpi/src/main.cpp $(MPI_OBJ_HYBRID) $(UTILS_OBJ_HYBRID) \
		$(LDFLAGS) $(LIBS) -o $@

run-mpi-hybrid: mpi_hybrid_pwr8
	OMP_NUM_THREADS=$(OMP_THREADS) mpirun -np $(NP) ./$(BUILD_DIR)/mpi_solver_hybrid $(ARGS)

mpi_hybrid_pwr8_openmpi: $(BUILD_DIR)/mpi_solver_hybrid_openmpi

$(BUILD_DIR)/mpi_solver_hybrid_openmpi: mpi/src/main.cpp $(MPI_OBJ_OPENMPI) $(UTILS_OBJ_OPENMPI) | $(BUILD_DIR)
	$(MPICXX_OPENMPI) $(CXXFLAGS_OPENMPI_HYBRID) mpi/src/main.cpp $(MPI_OBJ_OPENMPI) $(UTILS_OBJ_OPENMPI) \
		$(LDFLAGS) $(LIBS) -o $@

run-mpi-hybrid-openmpi: mpi_hybrid_pwr8_openmpi
	OMP_NUM_THREADS=$(OMP_THREADS) mpirun -np $(NP) ./$(BUILD_DIR)/mpi_solver_hybrid_openmpi $(ARGS)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -f $(UTILS_OBJ) $(SEQUENTIAL_OBJ) $(UTILS_OBJ_OMP) $(SEQUENTIAL_OBJ_OMP)
	rm -f $(MPI_OBJ) $(UTILS_OBJ_MPI)
	rm -f $(MPI_OBJ_PWR8) $(UTILS_OBJ_PWR8)
	rm -f $(MPI_OBJ_HYBRID) $(UTILS_OBJ_HYBRID)
	rm -f $(MPI_OBJ_OPENMPI) $(UTILS_OBJ_OPENMPI)
	rm -rf $(BUILD_DIR)
