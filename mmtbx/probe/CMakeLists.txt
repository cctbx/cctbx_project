# SCons is the build system that is used to build CCTBX projects.
# This CMakeLists.txt and the associated files in cmake/ were used before CCTBX native
# integration and during development to enable running the C++ code within the
# Visual Studio IDE for rapid debugging.

cmake_minimum_required(VERSION 3.10.0)
project(probe)

if (WIN32)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif ()

# Find the packages we require and can use
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(CCTBX REQUIRED)
find_package(Boost REQUIRED)

set(USE_PYTHON TRUE)
FIND_PACKAGE(PythonInterp)
if (PYTHONINTERP_FOUND)
  if (UNIX AND NOT APPLE)
    if (PYTHON_VERSION_MAJOR EQUAL 3)
        FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_SUFFIX})
        FIND_PACKAGE(PythonInterp 3)
        FIND_PACKAGE(PythonLibs 3 REQUIRED)
    else()
        FIND_PACKAGE(Boost COMPONENTS python)
        FIND_PACKAGE(PythonInterp)
        FIND_PACKAGE(PythonLibs REQUIRED)
    endif()
  else()
    if (PYTHON_VERSION_MAJOR EQUAL 3)
        FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
        FIND_PACKAGE(PythonInterp 3)
        FIND_PACKAGE(PythonLibs 3 REQUIRED)
    else()
        FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
        FIND_PACKAGE(PythonInterp)
        FIND_PACKAGE(PythonLibs REQUIRED)
    endif()
  endif()
else()
    set(USE_PYTHON FALSE)
    message("Python not found")
endif()
if (NOT Boost_FOUND)
    set(USE_PYTHON FALSE)
endif()

set(probe_SOURCES
  Scoring.cpp
  DotSpheres.cpp
  SpatialQuery.cpp
)

set(probe_HEADERS
  Common.h
  Scoring.h
  DotSpheres.h
  SpatialQuery.h
)

add_library(probeLib ${probe_SOURCES} ${probe_HEADERS})
target_include_directories(probeLib PUBLIC SYSTEM
  ${CCTBX_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIR})
target_link_libraries(probeLib
    PUBLIC
    ${CCTBX_LIBRARIES}
    ${Boost_LIBRARIES}
    ${PYTHON_LIBRARIES}
)
if (NOT WIN32)
  target_link_libraries(probeLib PUBLIC m)
endif()

add_executable(tst_probe tst_probe.cpp)
target_link_libraries(tst_probe probeLib)
if (NOT WIN32)
  target_link_libraries(tst_probe m)
endif()

if (USE_PYTHON)
  message(STATUS "Building Boost.Python library")
  message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
  message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
  message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
  message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")
  include_directories(../libpdb ../toolclasses ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})
  # Note: To make this work on Windows, you need to copy the boost_python*.dll file into the
  # executable directory to avoid an error about "DLL load failed: The specified module could not be found".
  # We also need to copy the test program into the Release (or other build) directory before running it.
  PYTHON_ADD_MODULE(mmtbx_probe_ext boost_python/probe_bpl.cpp)
  target_link_libraries(mmtbx_probe_ext PRIVATE probeLib)

  add_custom_command(
    TARGET mmtbx_probe_ext POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
      ${CMAKE_SOURCE_DIR}/tst_mmtbx_probe_ext.py
      ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()
