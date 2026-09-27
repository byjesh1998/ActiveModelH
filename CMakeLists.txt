cmake_minimum_required(VERSION 3.5.1)
project(activeH_heun)

set(CMAKE_CXX_STANDARD 14)

#set(CMAKE_CXX_FLAGS "-O3")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(Armadillo REQUIRED PATHS "~/armadillo")
find_package(Boost REQUIRED)
find_package(MKL REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

add_executable(activeH_heun main.cpp cfns.c)
target_link_libraries(activeH_heun mkl_rt ${ARMADILLO_LIBRARIES})
