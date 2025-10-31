CXX       := clang++
INCLUDES  := -Iutils/include -Isequential/include -Iparallel_openmp/include
LIBS      :=
LDFLAGS   :=

UTILS_SRC := utils/src/arg_parser.cpp utils/src/error_metrics.cpp
UTILS_OBJ := $(UTILS_SRC:.cpp=.o)

SEQUENTIAL_SRC := sequential/src/main.cpp sequential/src/solver.cpp
SEQUENTIAL_OBJ := $(filter-out sequential/src/main.o,$(SEQUENTIAL_SRC:.cpp=.o))

BUILD_DIR := build

CXXFLAGS_COMMON := -std=c++17 -Wall -Wextra $(INCLUDES)

CXXFLAGS_RELEASE := -O3 -march=native
CXXFLAGS_PROFILE := -O3 -g -gline-tables-only -fno-omit-frame-pointer -march=native

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


.PHONY: all clean test run-test sequential profile profile-deep debug trace run

all: sequential

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

GTEST_DIR     := /usr/local
GTEST_INCLUDE := -I$(GTEST_DIR)/include
GTEST_LIBS    := -L$(GTEST_DIR)/lib -lgtest -lgtest_main -pthread

test: test_runner

test_runner: $(UTILS_OBJ) utils/tests/test_arg_parser.cpp
	$(CXX) $(CXXFLAGS) $(GTEST_INCLUDE) \
		utils/tests/test_arg_parser.cpp $(UTILS_OBJ) \
		$(GTEST_LIBS) -o $@

run-test: test_runner
	./test_runner

sequential: $(BUILD_DIR)/sequential

$(BUILD_DIR)/sequential: sequential/src/main.cpp $(SEQUENTIAL_OBJ) $(UTILS_OBJ) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) sequential/src/main.cpp $(SEQUENTIAL_OBJ) $(UTILS_OBJ) \
		$(LDFLAGS) $(LIBS) -o $@
	$(POSTLINK)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

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

clean:
	rm -f $(UTILS_OBJ) $(SEQUENTIAL_OBJ) test_runner
	rm -rf $(BUILD_DIR)
