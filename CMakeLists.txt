## CMake build file for the chebfun library.
##
## This will build the library and its test suite. Compiled binaries will be
## written to build/${CMAKE_BUILD_TYPE}.
##
cmake_minimum_required(VERSION 2.8.12)

project(chebfun)

# Search for FFTW3
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(FFTW)
include_directories(${FFTW_INCLUDES})

set(LIB_FILES
    ${CMAKE_SOURCE_DIR}/src/chebfun.cpp
)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-std=c++11 -Wall")
#include_directories(src)
add_library(${PROJECT_NAME} ${LIB_FILES})

target_link_libraries(${PROJECT_NAME} PRIVATE ${FFTW_LIB})

function(install)
    # Disable install for test directory.
endfunction()
add_subdirectory(test)
