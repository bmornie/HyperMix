cmake_minimum_required(VERSION 3.18)

project(HyperMix LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)

execute_process(
    COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
    OUTPUT_VARIABLE pybind11_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE PYBIND11_RESULT
)

if(NOT PYBIND11_RESULT EQUAL 0)
    message(FATAL_ERROR "pybind11 not found. Install it using pip.")
endif()

find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(hypermix
    src/bindings.cpp
    src/hypergraph.cpp
    src/generator.cpp
    src/configuration_model.cpp
    utilities/hypergraph_motif.cpp
)

set_target_properties(hypermix PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/python"
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/python"

    LIBRARY_OUTPUT_DIRECTORY_RELEASE "${PROJECT_SOURCE_DIR}/python"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${PROJECT_SOURCE_DIR}/python"
)

target_include_directories(hypermix PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/utilities
)

if(MSVC)
    target_compile_options(hypermix PRIVATE /O2 /EHsc)
else()
    target_compile_options(hypermix PRIVATE -O3)
endif()