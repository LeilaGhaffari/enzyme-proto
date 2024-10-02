ENZYME_LIB ?=
CC = clang-18
CFLAGS = $(OPT) -Xclang -load -Xclang $(ENZYME_LIB) -Wall -Wextra -Wunused-variable -Wunused-function
LDFLAGS = -lm

CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -Wunused-variable -Wunused-function

# Source files
src.c = enzyme.c
src.cpp = adolc.cpp

# Object files
BUILD_DIR = build
obj.c = $(addprefix $(BUILD_DIR)/, $(src.c:.c=.o))
obj.cpp = $(addprefix $(BUILD_DIR)/, $(src.cpp:.cpp=.o))

# Executable names
TARGET_C = $(BUILD_DIR)/enzyme-exec
TARGET_CPP = $(BUILD_DIR)/adolc-exec

# Default rule
all: $(TARGET_C) $(TARGET_CPP)

# Link object files to create the C executable
$(TARGET_C): $(obj.c)
	$(CC) $(CFLAGS) -o $@ $(obj.c) $(LDFLAGS)

# Compile C source files to object files
$(BUILD_DIR)/%.o: %.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Link object files to create the C++ executable
$(TARGET_CPP): $(obj.cpp)
	$(CXX) $(CXXFLAGS) -I$(ADOLC_INCLUDE) -L$(ADOLC_LIB) -o $@ $(obj.cpp) -ladolc

# Compile C++ source files to object files
$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -I$(ADOLC_INCLUDE) -c $< -o $@

# Create the build directory if it doesn't exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Clean rule
clean:
	rm -rf $(BUILD_DIR)

.PHONY: all clean
