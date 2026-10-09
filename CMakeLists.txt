cmake_minimum_required(VERSION 3.16)
project(asap3sum LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_library(asap3sum
  src/thin_product.cpp
  src/exact_triangle.cpp
  src/three_sum.cpp)
target_include_directories(asap3sum PUBLIC include)
if(MSVC)
  target_compile_options(asap3sum PRIVATE /W4)
else()
  target_compile_options(asap3sum PRIVATE -Wall -Wextra -Wpedantic -Wshadow)
endif()

add_executable(asap3sum_demo examples/demo.cpp)
target_link_libraries(asap3sum_demo PRIVATE asap3sum)

enable_testing()
add_executable(asap3sum_tests tests/tests.cpp)
target_link_libraries(asap3sum_tests PRIVATE asap3sum)
add_test(NAME asap3sum_tests COMMAND asap3sum_tests)
