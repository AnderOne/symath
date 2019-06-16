cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 11)

project(SyMath++)

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

file(GLOB SYMATH_SRC "./src/symath/*.cpp")

add_library(
symath ${SYMATH_SRC}
)

add_executable(
demo ./src/test.cpp
)

target_include_directories(
symath PRIVATE
./inc
)

target_include_directories(
demo PRIVATE
./inc
)

target_link_libraries(
symath
"${GMPXX_LIB}"
"${GMP_LIB}"
)

target_link_libraries(
demo
symath
)