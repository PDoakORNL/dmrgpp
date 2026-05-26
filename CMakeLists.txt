cmake_minimum_required(VERSION 3.31)

project(
  dmrgpp
  VERSION 1.1
  LANGUAGES CXX)

option(BUILD_TESTING "Enable building tests" ON)
if(BUILD_TESTING)
  enable_testing()

  include(FetchContent)
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.11.0
    SYSTEM FIND_PACKAGE_ARGS NAMES Catch2)
  FetchContent_MakeAvailable(Catch2)
  include(Catch)
endif()

find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND)
  message(STATUS "Python3 found: ${Python3_EXECUTABLE} (${Python3_VERSION})")
else()
  message(STATUS "Python3 not found -- Python-based tests will be skipped")
endif()

include(GNUInstallDirs)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
if(APPLE)
  set(CMAKE_INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
else()
  set(CMAKE_INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()

add_subdirectory(util)
set(CHECK_VALUE_TARGET util::check_values)

add_subdirectory(PsimagLite)
add_subdirectory(dmrg)

add_subdirectory(cincuenta/src)
add_subdirectory(LanczosPlusPlus/src)
