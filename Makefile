CXX = clang++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2
INCLUDES = -Iutils/include -Isequential/include -Iparallel_openmp/include

BUILD_DIR = build

UTILS_SRC = utils/src/arg_parser.cpp utils/src/error_metrics.cpp
UTILS_OBJ = $(UTILS_SRC:.cpp=.o)

SEQUENTIAL_SRC = sequential/src/main.cpp sequential/src/solver.cpp
SEQUENTIAL_OBJ = $(filter-out sequential/src/main.o,$(SEQUENTIAL_SRC:.cpp=.o))

GTEST_DIR = /usr/local
GTEST_INCLUDE = -I$(GTEST_DIR)/include
GTEST_LIBS = -L$(GTEST_DIR)/lib -lgtest -lgtest_main -pthread

.PHONY: all clean test run-test

all:

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

test: test_runner

test_runner: $(UTILS_OBJ) utils/tests/test_arg_parser.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GTEST_INCLUDE) \
		utils/tests/test_arg_parser.cpp $(UTILS_OBJ) \
		$(GTEST_LIBS) -o test_runner

run-test: test_runner
	./test_runner

sequential: $(BUILD_DIR)/sequential

$(BUILD_DIR)/sequential: sequential/src/main.cpp $(SEQUENTIAL_OBJ) $(UTILS_OBJ) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) sequential/src/main.cpp $(SEQUENTIAL_OBJ) $(UTILS_OBJ) -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -f $(UTILS_OBJ)
	rm -f $(SEQUENTIAL_OBJ)
	rm -f test_runner
	rm -f sequential/src/*.o parallel_openmp/src/*.o
	rm -rf $(BUILD_DIR)
