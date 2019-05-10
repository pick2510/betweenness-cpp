cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

# Call cmake with -D TESTS=ON to set this flag to true.

project(DEM-Utils CXX C)

find_package(MPI REQUIRED)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()



set(Boost_USE_STATIC_LIBS        ON) # only find static libs
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
set(Boost_NO_BOOST_CMAKE    ON)
find_package(Boost REQUIRED COMPONENTS log_setup log serialization mpi filesystem)
include_directories(${Boost_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_PATH} ${PROJECT_SOURCE_DIR}/libeigen ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/inih )

set (CMAKE_CXX_STANDARD 14)
ADD_LIBRARY(sqlite3 STATIC
             ./src/sqlite3.c )


# Core and main are split. This allows us to link core to main and tests.

# Core library. *.cpp should be added here.

add_library(demtools STATIC ./src/utils.cpp ./src/graph.cpp ./src/graph_sqlite.cpp ./src/dumpfile.cpp ./src/decomposition.cpp src/potential_energy.cpp include/potential_energy.h)



# Main entry point.
add_executable(betweenness
  ./src/betweenness.cpp)


add_executable(grid
  ./src/grid.cpp)

add_executable(regrid
  ./src/regrid.cpp)

add_executable(betweenness_sqlite
  ./src/betweenness_sqlite.cpp)

add_executable(properties ./src/properties.cpp)

target_link_libraries(betweenness demtools ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} sqlite3)
target_link_libraries(regrid demtools ${Boost_LIBRARIES} sqlite3 dl pthread)
target_link_libraries(grid demtools ${Boost_LIBRARIES} sqlite3 dl pthread)
target_link_libraries(betweenness_sqlite demtools ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} sqlite3 dl pthread)
target_link_libraries(properties demtools ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} sqlite3 dl pthread)


# Add flags.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g")
set(CMAKE_C_FLAGS_DEBUG "-Og -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mtune=native")
set(CMAKE_C_FLAGS_RELEASE "-O3 -mtune=native")


#set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=undefined,leak")
#set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=undefined,leak")

#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg -no-pie")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg -no-pie")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
#SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg -no-pie")
