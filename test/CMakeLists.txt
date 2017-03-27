## CMake build file for the chebfun library test suite.
##
## This will build a Google Test test runner for all contained tests.
##
project(test_chebfun)

add_subdirectory(lib/gtest)

include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(${PROJECT_NAME} test_basic.cpp)
target_link_libraries(${PROJECT_NAME} chebfun gtest gtest_main)
