ENZYME_LIB=/home/leila/Enzyme/enzyme/build18g/Enzyme/ClangEnzyme-18.so
CC = clang-18
CFLAGS = $(OPT) -Xclang -load -Xclang $(ENZYME_LIB) -Wall -Wextra -Wunused-variable -Wunused-function
LDFLAGS = -lm

# Source files
SRCS = enzyme.c

# Object files
OBJS = $(SRCS:.c=.o)

# Executable name
TARGET = enzyme-exec

# Default rule
all: $(TARGET)

# Link object files to create the executable
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

# Compile source files to object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
