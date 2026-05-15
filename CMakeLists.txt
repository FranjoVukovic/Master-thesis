cmake_minimum_required(VERSION 3.15)

project(MasterThesis VERSION 1.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include(FetchContent)

# --- htslib ---
FetchContent_Declare(
    htslib
    GIT_REPOSITORY https://github.com/samtools/htslib.git
    GIT_TAG 1.20
)
FetchContent_MakeAvailable(htslib)

# --- yak ---
# yak uses a Makefile, not CMake, so we populate only and compile sources manually.
FetchContent_Declare(
    yak
    GIT_REPOSITORY https://github.com/lh3/yak.git
    GIT_TAG        v0.1
)
if(POLICY CMP0169)
    cmake_policy(SET CMP0169 OLD)
endif()
FetchContent_Populate(yak)

# --- thread_pool ---
find_package(thread_pool 4.0.0 QUIET)
if (NOT thread_pool_FOUND)
  FetchContent_Declare(
    thread_pool
    GIT_REPOSITORY https://github.com/rvaser/thread_pool
    GIT_TAG        4.0.0)
  FetchContent_GetProperties(thread_pool)
  if (NOT thread_pool_POPULATED)
    FetchContent_Populate(thread_pool)
    add_subdirectory(
      ${thread_pool_SOURCE_DIR}
      ${thread_pool_BINARY_DIR}
      EXCLUDE_FROM_ALL)
  endif ()
endif ()

file(GLOB YAK_SOURCES "${yak_SOURCE_DIR}/*.c")

# Exclude main.c to avoid multiple 'main' definitions
list(REMOVE_ITEM YAK_SOURCES "${yak_SOURCE_DIR}/main.c")

add_library(yak_lib STATIC ${YAK_SOURCES})

target_include_directories(yak_lib PUBLIC ${yak_SOURCE_DIR})

# yak's sources use pthreads and zlib
find_package(Threads REQUIRED)
target_link_libraries(yak_lib PUBLIC Threads::Threads z)

set(CORE_SOURCES
    src/Graph.cpp
    src/ContactMatrix.cpp
    src/BamParser.cpp
    src/TrioBinner.cpp
    src/BubbleDetector.cpp
    src/Phaser.cpp
    src/Chainer.cpp
    src/HamiltonianChainer.cpp
    src/HamiltonianPath.cpp
    src/UnzippedGraph.cpp
    src/ResultWriter.cpp
)

add_executable(MasterThesis src/main.cpp ${CORE_SOURCES})
target_include_directories(MasterThesis PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${htslib_SOURCE_DIR}
    ${yak_SOURCE_DIR}
)
target_link_libraries(MasterThesis PRIVATE hts yak_lib thread_pool::thread_pool)

enable_testing()

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.12.0
)
FetchContent_MakeAvailable(googletest)

add_executable(MasterThesis_Test test/test.cpp ${CORE_SOURCES})
target_include_directories(MasterThesis_Test PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${htslib_SOURCE_DIR}
    ${yak_SOURCE_DIR}
)

target_link_libraries(MasterThesis_Test PRIVATE hts yak_lib thread_pool::thread_pool gtest_main)

target_compile_definitions(MasterThesis_Test PRIVATE TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/test/data")
target_compile_options(MasterThesis_Test PRIVATE -fsanitize=address -fno-omit-frame-pointer -g)
target_link_options(MasterThesis_Test PRIVATE -fsanitize=address)

include(GoogleTest)
gtest_discover_tests(MasterThesis_Test)