## CMake build file for the chebfun library test suite.
##
## This will build a Google Test test runner for all contained tests.
##
project(test_chebfun)

# Search for FFTW3
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")
find_package(FFTW)
include_directories(${FFTW_INCLUDES})

add_subdirectory(lib/gtest)

include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(${PROJECT_NAME} test_basic.cpp)
target_link_libraries(${PROJECT_NAME} chebfun gtest gtest_main ${FFTW_LIB})
