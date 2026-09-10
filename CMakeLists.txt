cmake_minimum_required(VERSION 3.16)
project(geomath_tests LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(GEOMATH_REAL_TYPE "float" CACHE STRING "Real scalar type used by Geomath")
set_property(CACHE GEOMATH_REAL_TYPE PROPERTY STRINGS float double "long double")
if(NOT GEOMATH_REAL_TYPE STREQUAL "float" AND
   NOT GEOMATH_REAL_TYPE STREQUAL "double" AND
   NOT GEOMATH_REAL_TYPE STREQUAL "long double")
    message(FATAL_ERROR "GEOMATH_REAL_TYPE must be float, double, or long double")
endif()

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