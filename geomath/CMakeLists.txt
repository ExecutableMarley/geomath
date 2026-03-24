add_library(Geomath INTERFACE)

target_sources(Geomath INTERFACE
    Shapes/2D/Algorithms2D/Algorithms2D.cpp
    Shapes/2D/Algorithms2D/Delaunay/DelaunayTriangulation.cpp
    Shapes/3D/Algorithms3D/Algorithms3D.cpp
)

target_include_directories(Geomath INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(Geomath INTERFACE cxx_std_20)