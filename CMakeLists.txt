cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

# Call cmake with -D TESTS=ON to set this flag to true.
option(TESTS "build tests" OFF)

project(betweenness-cpp CXX C)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(Boost_USE_STATIC_LIBS        ON) # only find static libs
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
find_package(Boost COMPONENTS log)



# Core and main are split. This allows us to link core to main and tests.

# Core library. *.cpp should be added here.

# Main entry point.
add_executable(betweenness-cpp
  ./src/main.cpp ./src/utils.cpp ./src/graph.cpp)

# Link core to main.
target_link_libraries(betweenness-cpp ${Boost_LIBRARIES})
# Add flags.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O2")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

if(TESTS)

  enable_testing()

  # Process the CMakeLists.txt in third-party/googletest, sets project(gtest) and
  # project(gmock), gives access to ${gtest_SOURCE_DIR} and ${gmock_SOURCE_DIR}.
  add_subdirectory(third-party/googletest)

  include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
  include_directories(${gmock_SOURCE_DIR}/include ${gmock_SOURCE_DIR})

  # Tests. *-test.cpp should be added here.
  add_executable(tests
    ./test/main-test.cpp
    ./test/function-test.cpp
    ./test/class-test.cpp)

  # Link core, pthread and gtest to tests.
  target_link_libraries(tests
    pthread
    gtest
    gtest_main
    core)

endif()
