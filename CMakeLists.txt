cmake_minimum_required(VERSION 3.18)
project(AltAnn VERSION 0.3.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
option(ALTANN_RELEASE_BUILD "Bundle compiler runtimes for release packaging" OFF)

find_package(OpenMP QUIET COMPONENTS CXX)
if(NOT OpenMP_CXX_FOUND AND APPLE)
    # Apple Clang needs the separately installed OpenMP runtime. Support both
    # Homebrew prefixes, and allow a custom prefix through CMAKE_PREFIX_PATH.
    find_path(ALTANN_OMP_INCLUDE omp.h
        HINTS /opt/homebrew/opt/libomp/include /usr/local/opt/libomp/include)
    find_library(ALTANN_OMP_LIBRARY omp
        HINTS /opt/homebrew/opt/libomp/lib /usr/local/opt/libomp/lib)
    if(ALTANN_OMP_INCLUDE AND ALTANN_OMP_LIBRARY)
        add_library(OpenMP::OpenMP_CXX INTERFACE IMPORTED)
        set_target_properties(OpenMP::OpenMP_CXX PROPERTIES
            INTERFACE_COMPILE_OPTIONS "-Xpreprocessor;-fopenmp"
            INTERFACE_INCLUDE_DIRECTORIES "${ALTANN_OMP_INCLUDE}"
            INTERFACE_LINK_LIBRARIES "${ALTANN_OMP_LIBRARY}")
    endif()
endif()
if(NOT TARGET OpenMP::OpenMP_CXX)
    message(FATAL_ERROR "OpenMP is required. On macOS, install libomp with Homebrew.")
endif()

add_executable(altann-core
    src/main.cpp src/loaders.cpp src/model.cpp src/global.cpp
    src/local.cpp src/output.cpp src/output_files.cpp)
target_link_libraries(altann-core PRIVATE OpenMP::OpenMP_CXX)
target_compile_options(altann-core PRIVATE -Wall -Wextra)
# Preserve the runtime search path for a user-installed OpenMP library.
if(APPLE)
    set_target_properties(altann-core PROPERTIES INSTALL_RPATH "@loader_path/../lib"
        INSTALL_RPATH_USE_LINK_PATH TRUE)
else()
    set_target_properties(altann-core PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")
    if(ALTANN_RELEASE_BUILD AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_link_options(altann-core PRIVATE -static-libstdc++ -static-libgcc)
    endif()
endif()

install(TARGETS altann-core RUNTIME DESTINATION bin)
install(PROGRAMS bin/altann DESTINATION bin)
install(DIRECTORY altann/ DESTINATION share/altann/altann FILES_MATCHING
    PATTERN "*.py" PATTERN "__pycache__" EXCLUDE)
install(FILES LICENSE THIRD_PARTY.md AUTHORS.md DESTINATION share/altann)

include(CTest)
if(BUILD_TESTING)
    find_package(Python3 3.9 REQUIRED COMPONENTS Interpreter)
    add_test(NAME integration COMMAND ${CMAKE_COMMAND} -E env
        "ALTANN_CORE=$<TARGET_FILE:altann-core>"
        ${Python3_EXECUTABLE} -m unittest discover -s tests -v)
    set_tests_properties(integration PROPERTIES WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
endif()
