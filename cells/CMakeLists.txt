# -----------------------------------------------------------------------------
#
# Copyright (C) 2021 CERN & University of Surrey for the benefit of the
# BioDynaMo collaboration. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
#
# See the LICENSE file distributed with this work for details.
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# -----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.19.3)
project(cells)

# BioDynaMo currently uses the C++17 standard.
set(CMAKE_CXX_STANDARD 17)

# Use BioDynaMo in this project.
find_package(BioDynaMo REQUIRED)

# Find preCICE
find_package(precice REQUIRED)

# See UseBioDynaMo.cmake in your BioDynaMo build folder for details.
include(${BDM_USE_FILE})

# Consider all files in src/ for BioDynaMo simulation.
include_directories("src")
include_directories(${precice_INCLUDE_DIRS})  # Ensure preCICE headers are found

file(GLOB_RECURSE PROJECT_HEADERS src/*.h)
file(GLOB_RECURSE PROJECT_SOURCES src/*.cc)

# Add preCICE library directly to bdm_add_executable
bdm_add_executable(${CMAKE_PROJECT_NAME}
                   HEADERS ${PROJECT_HEADERS}
                   SOURCES ${PROJECT_SOURCES}
                   LIBRARIES ${BDM_REQUIRED_LIBRARIES} precice::precice)

# Consider all files in test/ for GoogleTests.
include_directories("test")
file(GLOB_RECURSE TEST_SOURCES test/*.cc)
file(GLOB_RECURSE TEST_HEADERS test/*.h)

bdm_add_test(${CMAKE_PROJECT_NAME}-test
             SOURCES ${TEST_SOURCES}
             HEADERS ${TEST_HEADERS}
             LIBRARIES ${BDM_REQUIRED_LIBRARIES} ${CMAKE_PROJECT_NAME})