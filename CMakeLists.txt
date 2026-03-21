cmake_minimum_required(VERSION 3.15)

project(MasterThesis VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Outputs the executable to a 'bin' folder in your build directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Assuming main.cpp is in the root directory
add_executable(MasterThesis src/main.cpp)