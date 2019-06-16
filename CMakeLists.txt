cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 11)

project(SyMath++)

add_library(
symath
func_tree_cpy.cpp
func_tree_dif.cpp
func_tree_get.cpp
func_tree_red.cpp
func_tree_str.cpp
func_tree.cpp
long.cpp
)

add_executable(
demo test.cpp
)

find_library(GMPXX_LIB libgmpxx.a)
find_library(GMP_LIB libgmp.a)

find_file(GMP_INC gmpxx.h)

if (NOT GMPXX_LIB)
message(FATAL_ERROR "File 'libgmpxx.a' not found! GMP library must be installed!")
endif()

if (NOT GMP_LIB)
message(FATAL_ERROR "File 'libgmp.a' not found! GMP library must be installed!")
endif()

if (NOT GMP_INC)
message(FATAL_ERROR "File <gmpxx.h> not found! GMP library must be installed!")
endif()

target_link_libraries(
symath
"${GMPXX_LIB}"
"${GMP_LIB}"
)

target_link_libraries(
demo symath
)
