cmake_minimum_required(VERSION 3.5)
project(Aeroland++)

# C++14
set(CMAKE_CXX_STANDARD 14)

# Setup C++ Compile Flags
add_compile_options(-Wall -Wextra -Wpedantic -Wno-unused-parameter -Werror=return-type -O1)

find_package(Boost REQUIRED)   # for odeint
find_package(OpenMP REQUIRED)  # for odeint

# Bundle all the flags required to build OpenMP
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

# Build ECOS shared library
execute_process(COMMAND make shared WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/ecos")

# Include ECOS shared library
include_directories(ecos/include)
include_directories(ecos/external/SuiteSparse_config)

# For vector definitions
include_directories(SYSTEM eigen)

# Link ECOS library for use in problem.c
link_directories(ecos)
link_libraries(ecos)

# Target
add_executable(Aeroland++ src/main.cpp include/problem.h src/problem.c include/successive.h)
target_link_libraries(Aeroland++ m)  # for math.h
