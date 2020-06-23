# include custom add_test function stolen from here: 
# https://eb2.co/blog/2015/06/driving-boost.test-with-cmake/
include(BoostTestHelpers)

# find boost unit test framework
find_package(Boost COMPONENTS unit_test_framework REQUIRED)

# add boost include directory to look for headers
include_directories( ${Boost_INCLUDE_DIR} )

# copy file to the test directory
configure_file(normal_0_1_data.csv ../test/normal_0_1_data.csv COPYONLY)

# needed to use boost unit_test dynamic library to generate main function
# for tests. FIXME: Problematic upon compiler switches!
add_definitions(-DBOOST_TEST_DYN_LINK) 
add_boost_test(algorithms_test.cpp StatLab)
