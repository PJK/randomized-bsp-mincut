cmake_minimum_required (VERSION 2.8)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
project(parallel-minimum-cut)
include(CTest)

enable_testing()

set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /opt/local/lib)

set(CTEST_MEMORYCHECK_COMMAND "/usr/bin/valgrind")
set(MEMORYCHECK_COMMAND_OPTIONS "")

add_compile_options(-Wall)

# Ancient g++ on Euler needs this for the chrono namespace
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-deprecated-declarations")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -ggdb -g -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -march=native -flto -DNDEBUG")

if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -xavx -ipo")
endif()

find_package(MPI)

if (APPLE AND (NOT MPI_CXX_LIBRARIES OR NOT MPI_C_LIBRARIES)) 
	set (MPI_C_COMPILER mpicc-openmpi-mp)
	set (MPI_CXX_COMPILER mpicxx-openmpi-mp)
endif(APPLE AND (NOT MPI_CXX_LIBRARIES OR NOT MPI_C_LIBRARIES))

find_package(MPI REQUIRED)

set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
include_directories(${MPI_INCLUDE_PATH})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

option(PROFILE_TIMING "Print timing data" OFF)
if(PROFILE_TIMING)
	add_definitions(-DPROFILE_TIMING)
endif(PROFILE_TIMING)

option(DAINT_DEBUG "Print Daint crashes RC debugging info" OFF)
if(DAINT_DEBUG)
	add_definitions(-DDAINT_DEBUG)
endif(DAINT_DEBUG)

include_directories(src)
include_directories(lib)
subdirs(src)
