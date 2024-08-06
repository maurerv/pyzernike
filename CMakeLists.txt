cmake_minimum_required(VERSION 3.15...3.27)

project(
    pyzernike
    VERSION ${SKBUILD_PROJECT_VERSION}
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm|i386")
    set(ARCH_FLAGS "-O3 -ftree-vectorize")
else()
    set(ARCH_FLAGS "-O3 -march=native")
endif()
set(COMMON_FLAGS "-funroll-loops -ffast-math -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ARCH_FLAGS} ${COMMON_FLAGS}")


find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)
pybind11_add_module(_core src/bindings.cpp)
target_include_directories(_core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/lib)

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Zernike)
#pybind11_add_module(_core src/bindings.cpp)
#target_include_directories(_core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Zernike)
#target_link_libraries(_core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Zernike)

install(TARGETS _core DESTINATION pyzernike)