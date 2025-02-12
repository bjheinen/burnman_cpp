# Target Executable (to make)
TARGET_EXEC := burnman

# Compiler
CXX := g++

# Directories
CWD := $(shell pwd)
SRC_DIR := $(CWD)/src
INCLUDE_DIR := $(CWD)/include
BUILD_DIR := $(CWD)/build
BIN_DIR := $(CWD)/bin

# For GSL, etc. (Change if non-standard)
# Only defined here so they can be changed
# both should be on the system path anyway
EXTRA_INCLUDE := /usr/local/include
EXTRA_LIB := /usr/local/lib
LDFLAGS := -lgsl -lgslcblas -lm

# Expand to paths
INCLUDE_FLAGS := -I$(INCLUDE_DIR) $(addprefix -I, $(EXTRA_INCLUDE))
LIB_PATHS := $(addprefix -L, $(EXTRA_LIB))

# Compiler flags
CXXFLAGS := -Wall -Wextra -std=c++17
CPPFLAGS := -MMD -MP $(INCLUDE_FLAGS)

# Find source files
SOURCES := $(wildcard $(SRC_DIR)/**/*.cpp)
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SOURCES))
DEPENDS := $(OBJECTS:.o=.d)

# Output target
TARGET := $(BIN_DIR)/$(TARGET_EXEC)

# Default rule
all: $(TARGET)

# Link and build
$(TARGET) : $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(LIB_PATHS) $^ -o $@ $(LDFLAGS)

# Compile source files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# Create build dirs
$(BUILD_DIR):
	mkdir -p $@

$(BIN_DIR):
	mkdir -p $@

-include $(DEPENDS)

# Clean up
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

# Maybe debug/release with opt flags?
# Debug and Release builds
#debug: CXXFLAGS += -g -O0
#debug: $(TARGET)

#release: CXXFLAGS += -O2
#release: $(TARGET)

.PHONY: all clean
	
