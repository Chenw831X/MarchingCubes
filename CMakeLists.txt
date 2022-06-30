cmake_minimum_required(VERSION 3.22)
project(MC)

option(BUILD_TESTING "Build tests" ON)

set(CMAKE_CXX_STANDARD 17)
# set compile flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -march=native")

# project source files
file(GLOB SRCFILES src/*.cpp)

add_library(${PROJECT_NAME}_dev ${SRCFILES})

target_include_directories(${PROJECT_NAME}_dev PUBLIC src)

# Eigen, igl, spdlog
target_include_directories(${PROJECT_NAME}_dev PUBLIC /usr/local/include)
# disable Eigen's multi-threading in multi-threaded application
target_compile_definitions(${PROJECT_NAME}_dev PUBLIC EIGEN_DONT_PARALLELIZE)

# add -pthread to compilation and linking
find_package(Threads REQUIRED)
target_link_libraries(${PROJECT_NAME}_dev PUBLIC ${CMAKE_THREAD_LIBS_INIT})

find_package(Boost COMPONENTS filesystem REQUIRED)
target_link_libraries(${PROJECT_NAME}_dev PUBLIC ${Boost_LIBRARIES})


# Tests
if(BUILD_TESTING)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests/sciVis)
endif()

enable_testing()
if(BUILD_TESTING)
    add_test(NAME test_sciVis COMMAND $<TARGET_FILE:test_sciVis>)
endif()
