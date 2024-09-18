# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++11


# Source files
SRCS = src/main.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Executable name
TARGET = exec

# Default rule
all: $(TARGET)

# Link object files to create the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -I$(ADOLC_INCLUDE) -Iinclude -L$(ADOLC_LIB) -o $(TARGET) $(OBJS) -ladolc

# Compile source files to object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I$(ADOLC_INCLUDE) -Iinclude -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
