# include custom add_test function stolen from here: 
# https://eb2.co/blog/2015/06/driving-boost.test-with-cmake/
include(BoostTestHelpers)

# find boost unit test framework
find_package(Boost COMPONENTS unit_test_framework REQUIRED)

# add boost include directory to look for headers
include_directories( ${Boost_INCLUDE_DIR} )

set(CMAKE_CXX_STANDARD 20)

# copy files to test dir 
configure_file(normal.csv normal.csv COPYONLY)
configure_file(uniform.csv uniform.csv COPYONLY)

# add executable for the test
add_executable(statistics_test statistics_test.cpp)

# needed to use boost unit_test dynamic library to generate main function
# for tests.
# add_definitions(-DBOOST_TEST_DYN_LINK) 
add_test(statistics_test statistics_test)


