# Target Executable (to make)
TARGET_EXEC := burnman
TEST_EXEC := run_tests

# Compiler
CXX := g++

# Directories
CWD := $(shell pwd)
SRC_DIR := $(CWD)/src
INCLUDE_DIR := $(CWD)/include
BUILD_DIR := $(CWD)/build
BIN_DIR := $(CWD)/bin
TEST_DIR := $(CWD)/tests
INCLUDE_DIR_TEST := $(TEST_DIR)/include

# For GSL, etc. (Change if non-standard)
# Only defined here so they can be changed
# both should be on the system path anyway
EXTRA_INCLUDE := /usr/local/include
EXTRA_LIB := /usr/local/lib
LDFLAGS_COMMON := -lgsl -lgslcblas -lm
LDFLAGS_TEST := -lCatch2Main -lCatch2

# Expand to paths
INCLUDE_FLAGS_COMMON := -I$(INCLUDE_DIR) $(addprefix -I, $(EXTRA_INCLUDE))
INCLUDE_FLAGS_TEST := -I$(INCLUDE_DIR_TEST)
LIB_PATHS := $(addprefix -L, $(EXTRA_LIB))

# Compiler flags
CXXFLAGS := -Wall -Wextra -pedantic -Wshadow -Wconversion -std=c++17
CPPFLAGS := -MMD -MP $(INCLUDE_FLAGS_COMMON)

# Build type specific flags
DEBUG_FLAGS = -g -Og
RELEASE_FLAGS = -O3 -DNDEBUG
TEST_FLAGS = -g -O1 -fno-omit-frame-pointer #-fsanitize=address,undefined

# Select Mode (default: release)
BUILD_MODE ?= release
ifeq ($(BUILD_MODE), debug)
  CXXFLAGS += $(DEBUG_FLAGS)
  BUILD_SUFFIX = _debug
else ifeq ($(BUILD_MODE), test)
  CXXFLAGS += $(TEST_FLAGS)
  BUILD_SUFFIX = _test
else
  CXXFLAGS += $(RELEASE_FLAGS)
  BIN_SUFFIX =
endif

# Find source files
SOURCES := $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(SRC_DIR)/**/*.cpp)
TEST_SOURCES := $(wildcard $(TEST_DIR)/*.cpp) $(wildcard $(TEST_DIR)/**/*.cpp)

# Generate object fnames
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SOURCES))
TEST_OBJECTS := $(patsubst $(TEST_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(TEST_SOURCES))

# Generate dependencies
DEPENDS := $(OBJECTS:.o=.d)

# Output target
TARGET := $(BIN_DIR)/$(TARGET_EXEC)$(BUILD_SUFFIX)
TEST_TARGET := $(BIN_DIR)/$(TEST_EXEC)$(BUILD_SUFFIX)

# Default rule
all: $(TARGET)

# Link and build
$(TARGET) : $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(LIB_PATHS) $^ -o $@ $(LDFLAGS_COMMON)

# Link and build test executable
$(TEST_TARGET): $(OBJECTS) $(TEST_OBJECTS) | $(BIN_DIR)
	$(CXX) $(LIB_PATHS) $^ -o $@ $(LDFLAGS_COMMON) $(LDFLAGS_TEST)

# Compile source files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# Compile tests
$(BUILD_DIR)/%.o: $(TEST_DIR)/%.cpp | $(BUILD_DIR)
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(INCLUDE_FLAGS_TEST) $(CXXFLAGS) -c $< -o $@

# Create build dirs
$(BUILD_DIR):
	mkdir -p $@

$(BIN_DIR):
	mkdir -p $@

-include $(DEPENDS)

# Run tests
test: $(TEST_TARGET)
	$(TEST_TARGET)

# Clean up
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

# Maybe debug/release with opt flags?
# Debug and Release builds
#debug: CXXFLAGS += -g -O0
#debug: $(TARGET)

#release: CXXFLAGS += -O2
#release: $(TARGET)

.PHONY: all clean test
	
