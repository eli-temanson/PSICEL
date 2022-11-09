# Check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION 3.0 FATAL_ERROR)

set(PROJECT_NAME ic_reso)
project(${PROJECT_NAME})
add_subdirectory(vendor/catima)
add_subdirectory(vendor/yaml-cpp)
add_subdirectory(vendor/spdlog)

find_package(ROOT)

set(CMAKE_CXX_FLAGS "${ROOT_CXX_FLAGS}")

include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})

link_directories(${ROOT_LIBRARY_DIR})

file(GLOB SOURCES "./src/*.cpp")
add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} ${ROOT_LIBRARIES} catima yaml-cpp spdlog)

set_target_properties(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${MY_BINARY_DIR})

target_compile_definitions(${PROJECT_NAME} PRIVATE YAML_CPP_STATIC_DEFINE)
