cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 20)

project(psicel)

set(MY_BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(MY_LIBRARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)

add_subdirectory(psicel)
