ENZYME_LIB ?=
ADOLC_INCLUDE ?=
ADOLC_LIB ?=

CC = clang-18
CFLAGS = $(OPT) -Wall -Wextra -Wunused-variable -Wunused-function -Iinclude
LDFLAGS = -lm

ifneq ($(ENZYME_LIB),)
    CFLAGS += -Xclang -load -Xclang $(ENZYME_LIB)
endif

CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -Wunused-variable -Wunused-function -Iinclude -I$(ADOLC_INCLUDE)

# Directories
SRCDIR = src
INCDIR = include
BUILDDIR = build

# Source files
SRC_C = $(SRCDIR)/enzyme.c
SRC_CPP = $(SRCDIR)/adolc.cpp

# Object files
OBJ_C = $(BUILDDIR)/$(notdir $(SRC_C:.c=.o))
OBJ_CPP = $(BUILDDIR)/$(notdir $(SRC_CPP:.cpp=.o))

# Executable names
TARGET_C = $(BUILDDIR)/enzyme-exec
TARGET_CPP = $(BUILDDIR)/adolc-exec

# Default rule
all: $(TARGET_C) $(TARGET_CPP)

# Link object files to create the executable (C)
$(TARGET_C): $(OBJ_C)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Link object files to create the executable (C++)
$(TARGET_CPP): $(OBJ_CPP)
	$(CXX) $(CXXFLAGS) -I$(ADOLC_INCLUDE) -I$(INCDIR) -L$(ADOLC_LIB) -o $@ $^ -ladolc

# Compile source files to object files (C)
$(BUILDDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

# Compile source files to object files (C++)
$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -I$(ADOLC_INCLUDE) -I$(INCDIR) -c $< -o $@

# Clean rule
clean:
	rm -f $(BUILDDIR)/*.o $(TARGET_C) $(TARGET_CPP)

.PHONY: all clean
