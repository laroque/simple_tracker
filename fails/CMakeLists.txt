cmake_minimum_required(VERSION 2.6.4)
set(CMAKE_C_COMPILER "/usr/bin/gcc")
set(CMAKE_CXX_COMPILER "/usr/bin/g++")
set(CMAKE_VERBOSE_MAKEFILE on)
project(CyclotronTracker)
add_executable(hello hello.cpp)

find_package(gsl)
INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ${gsl_INCLUDE_DIRS})
LINK_DIRECTORIES(${LINK_DIRECTORIES} ${gsl_LINK_DIRS})

add_executable(tracker tracker_aspen.cpp)
target_link_libraries(tracker gsl)
