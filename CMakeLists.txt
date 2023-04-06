# Copyright 2023.01.13 Jaewoo An in Korea Institute of sicience and technology
#
# It is for build C++ code with MuJoCo simulation
# You may not use this file without author`s permission
#

# CMake minimun version
cmake_minimum_required(VERSION 3.16)

# Project information
project(
  Simulation 
  VERSION 2.3.1
  DESCRIPTION "The basic MuJoCo simulate with Eigen" 
  LANGUAGES CXX C
)

# cmake .. -DCMAKE_BUILD_TYPE=Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)        # LTO option
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_CXX_FLAGS} -Wl,-no-as-needed -fexceptions")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")           # -Wall -Wextra -Wpedantic -Werror")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -march=native")   # -Ofast, -O3, -O2, -O1, -Os


find_package(X11 REQUIRED)
find_package(Eigen3 REQUIRED)
find_library(GLFW libglfw.so.3 REQUIRED) # HINTS ${CMAKE_SOURCE_DIR}/lib/glfw/glfw-3.3.8)

message(STATUS "PROJECT_NAME: " ${PROJECT_NAME})
message(STATUS "MUJOCO_HOME: " ${CMAKE_SOURCE_DIR})
message(STATUS "GLFW lib found at: " ${GLFW})

# Make exe file
add_executable(${PROJECT_NAME}
  ${CMAKE_SOURCE_DIR}/lib/mujoco-2.3.1/simulate/glfw_dispatch.cc
  ${CMAKE_SOURCE_DIR}/lib/mujoco-2.3.1/simulate/simulate.cc
  ${CMAKE_SOURCE_DIR}/lib/mujoco-2.3.1/simulate/uitools.cc
  ${CMAKE_SOURCE_DIR}/lib/lodepng/lodepng-master/lodepng.cpp
  ${CMAKE_SOURCE_DIR}/include/ARBMLlib/ARBML.cpp
  ${CMAKE_SOURCE_DIR}/include/ARBMLlib/RigidBody.cpp
  ${CMAKE_SOURCE_DIR}/src/Robot_Control.cpp
  ${CMAKE_SOURCE_DIR}/src/main.cpp
)

# include header files in target file
target_include_directories(${PROJECT_NAME} PUBLIC
  ${CMAKE_SOURCE_DIR}/lib/mujoco-2.3.1/simulate
  ${CMAKE_SOURCE_DIR}/lib/mujoco-2.3.1/include
  ${CMAKE_SOURCE_DIR}/lib/lodepng/lodepng-master
  ${CMAKE_SOURCE_DIR}/lib/glfw/glfw-3.3.8/include
  ${CMAKE_SOURCE_DIR}/include
  ${EIGEN3_INCLUDE_DIR}
)

target_link_libraries(${PROJECT_NAME} PUBLIC
  libmujoco.so.2.3.1
  -lpthread
  # -lglfw3
  # lodepng
  ${CMAKE_DL_LIBS}
  ${GLFW}
  ${X11_LIBRARIES}
)

target_link_directories(${PROJECT_NAME} PUBLIC
  lib/mujoco-2.3.1/lib
)