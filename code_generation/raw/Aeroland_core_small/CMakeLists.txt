cmake_minimum_required(VERSION 3.5)
project(Aeroland_core_small C)

set(CMAKE_C_STANDARD 99)

# Build ECOS shared library
execute_process(COMMAND make shared WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/ecos")

include_directories(ecos/include)
include_directories(ecos/external/SuiteSparse_config)

link_directories(ecos)
link_libraries(ecos)
