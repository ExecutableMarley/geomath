add_executable(geomath_tests)

file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS "test_*.cpp")

target_sources(geomath_tests PRIVATE
    main.cpp
    ${TEST_SOURCES}
)

# doctest forward-declares std::tuple in some builds. use the real STL headers on MSVC
# to avoid warning C5285/N5014 on the Windows runner.
target_compile_definitions(geomath_tests PRIVATE DOCTEST_CONFIG_USE_STD_HEADERS)

target_link_libraries(geomath_tests PRIVATE Geomath)

# Register doctest binary
add_test(
    NAME geomath_tests
    COMMAND geomath_tests
)