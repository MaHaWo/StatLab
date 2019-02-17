# --- TESTS ---
# All tests registered with 'add_test' will be run when calling 'make test'
#
# NOTE: CTest does *not* compile its test executables. Running 'make test'
# *without* compiling the appropriate executable will cause the test to fail!
#
# add test as regular executable
add_executable(statistics_test statistics_test.cc)
target_compile_options(statistics_test
    PRIVATE
        $<$<CONFIG:DEBUG>:--coverage>
)
# link the library to it
target_link_libraries(statistics_test
    # can stay private, this target is never exported
    PRIVATE
        statistics
        $<$<CONFIG:DEBUG>:--coverage>
)

# register the test in CTest
add_test(NAME test-statistics
    # automatically replaces target name with executable location
    COMMAND statistics_test
)
