cmake_minimum_required(VERSION 3.10)

project(SyMath++)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_CXX_STANDARD 11)

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

add_library(symath ${SYMATH_SRC})
target_include_directories(symath PRIVATE ./inc)
target_link_libraries(
symath "${GMPXX_LIB}" "${GMP_LIB}"
)

add_executable(demo ./src/demo.cpp)
target_include_directories(demo PRIVATE ./inc)
target_link_libraries(demo symath)

set(Boost_USE_STATIC_LIBS ON)

find_package(
Boost COMPONENTS unit_test_framework
)

if (Boost_FOUND)

add_executable(test ./src/test.cpp)
target_include_directories(
test PRIVATE ./inc
)
target_link_libraries(
test
${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
target_link_libraries(
test symath
)

endif()
