add_executable(geomath_tests)

file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS "test_*.cpp")

target_sources(geomath_tests PRIVATE
    main.cpp
    ${TEST_SOURCES}
)

target_link_libraries(geomath_tests PRIVATE Geomath)

# Register doctest binary
add_test(
    NAME geomath_tests
    COMMAND geomath_tests
)