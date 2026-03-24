cmake_minimum_required(VERSION 3.16)
project(geomath_tests LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_testing()

add_subdirectory(geomath)
add_subdirectory(tests)

# Amalgamation target
set(AMALGAMATE_EXECUTABLE
    ${CMAKE_SOURCE_DIR}/tools/amalgamate
)

set(DIST_DIR
    ${CMAKE_SOURCE_DIR}/dist
)

add_custom_target(amalgamate
    COMMAND ${CMAKE_COMMAND} -E make_directory ${DIST_DIR}

    COMMAND ${AMALGAMATE_EXECUTABLE}
        -i ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/Math.hpp
        ${DIST_DIR}/geomath.hpp

    COMMAND ${AMALGAMATE_EXECUTABLE}
        -i ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/MathAmalgamated.cpp
        ${DIST_DIR}/geomath.cpp

    COMMENT "Amalgamating geomath headers and sources"
)

# Combined target
add_custom_target(test_and_amalgamate
    COMMAND ${CMAKE_CTEST_COMMAND}
            --output-on-failure
    COMMAND ${CMAKE_COMMAND} --build . --target amalgamate

    DEPENDS geomath_tests
    COMMENT "Running tests and amalgamating on success"
)