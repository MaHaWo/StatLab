cmake_minimum_required(VERSION 3.10)

# add path where to look for modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeFiles/")

# add custom test command for verbose output
add_custom_target(test_all COMMAND ${CMAKE_CTEST_COMMAND}  --output-on-failure)

# set up header only lib
project(StatLab VERSION 0.1 DESCRIPTION "Online, stable  algorithms for simple statistics")
add_library(StatLab INTERFACE)
target_include_directories(StatLab INTERFACE include/)

# set c++17 standard
set (CMAKE_CXX_STANDARD 17)

# enable testing and shit
ENABLE_TESTING()
ADD_SUBDIRECTORY( test )
