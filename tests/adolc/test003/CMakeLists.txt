# Set CMake version
cmake_minimum_required(VERSION 3.14)

# Set project name
set(project_name "ADOLC_MINIMAL")
project(${project_name})

# Set path to header files directories
include_directories($ENV{ADOLC_INCLUDE})

# Set path to executable directories
link_directories($ENV{ADOLC_LIB})

# Add source files to compile to the project
add_executable(exec main.cpp)

# Add external libraries
target_link_libraries(exec -ladolc)
