cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

# Call cmake with -D TESTS=ON to set this flag to true.

project(DEM-Utils CXX C)

find_package(MPI REQUIRED)

#find_package(OpenMP)
#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/SQLiteCpp)


set(Boost_USE_STATIC_LIBS        ON) # only find static libs
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
set(Boost_NO_BOOST_CMAKE    ON)
find_package(Boost COMPONENTS log_setup log serialization mpi filesystem)
include_directories(${Boost_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_PATH} ${PROJECT_SOURCE_DIR}/libeigen ${PROJECT_SOURCE_DIR}/include  ${CMAKE_CURRENT_LIST_DIR}/SQLiteCpp)

set (CMAKE_CXX_STANDARD 11)



# Core and main are split. This allows us to link core to main and tests.

# Core library. *.cpp should be added here.

# Main entry point.
add_executable(betweenness-cpp
  ./src/main.cpp ./src/utils.cpp ./src/graph.cpp)




target_link_libraries(betweenness-cpp ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
# Add flags.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g")
set(CMAKE_C_FLAGS_DEBUG "-Og -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mtune=native")
set(CMAKE_C_FLAGS_RELEASE "-O3 -mtune=native")
