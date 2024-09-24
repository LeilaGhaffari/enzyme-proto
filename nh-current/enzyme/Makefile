ENZYME_LIB=/home/leila/Enzyme/enzyme/build18g/Enzyme/ClangEnzyme-18.so
CC = clang-18
CFLAGS = $(OPT) -Xclang -load -Xclang $(ENZYME_LIB) -Wall -Wextra -Wunused-variable -Wunused-function
LDFLAGS = -lm

# Source files
SRCS = src/main.c

# Object files
OBJS = $(SRCS:.c=.o)

# Executable name
TARGET = exec

# Default rule
all: $(TARGET)

# Link object files to create the executable
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -Iinclude -o $(TARGET) $(OBJS) $(LDFLAGS)

# Compile source files to object files
%.o: %.c
	$(CC) $(CFLAGS) -Iinclude -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
