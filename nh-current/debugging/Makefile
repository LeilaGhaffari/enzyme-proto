ENZYME_LIB=/home/leila/Enzyme/enzyme/build18g/Enzyme/ClangEnzyme-18.so
CC = clang-18
CFLAGS = $(OPT) -Xclang -load -Xclang $(ENZYME_LIB) -Wall -Wextra -Wunused-variable -Wunused-function
LDFLAGS = -lm

CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -Wunused-variable -Wunused-function

# Source files
src.c = enzyme.c
src.cpp = adolc.cpp

# Object files
obj.c = $(src.c:.c=.o)
obj.cpp = $(src.cpp:.cpp=.o)

# Executable names
TARGET_C = enzyme-exec
TARGET_CPP = adolc-exec

# Default rule
all: $(TARGET_C) $(TARGET_CPP)

# Link object files to create the C executable
$(TARGET_C): $(obj.c)
	$(CC) $(CFLAGS) -o $(TARGET_C) $(obj.c) $(LDFLAGS)

# Compile C source files to object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Link object files to create the C++ executable
$(TARGET_CPP): $(obj.cpp)
	$(CXX) $(CXXFLAGS) -I$(ADOLC_INCLUDE) -L$(ADOLC_LIB) -o $(TARGET_CPP) $(obj.cpp) -ladolc

# Compile C++ source files to object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I$(ADOLC_INCLUDE) -c $< -o $@

# Clean rule
clean:
	rm -f $(obj.c) $(obj.cpp) $(TARGET_C) $(TARGET_CPP)

.PHONY: all clean

