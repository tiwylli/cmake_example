cmake_minimum_required(VERSION 3.4...3.18)
project(diffvoronoi)

add_subdirectory(3rdparty/geogram)

add_subdirectory(pybind11)
pybind11_add_module(diffvoronoi src/main.cpp)

include_directories(3rdparty/geogram/src/lib)
target_link_libraries(diffvoronoi PRIVATE geogram)

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(diffvoronoi
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
