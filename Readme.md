# GeoMath

A lightweight C++ geometry and math library focused on shape
abstractions and interaction tests.

GeoMath is distributed as two amalgamated files and can be dropped
directly into projects with minimal setup.

## Features

## Integration

Simply add the amalgamated header & source file into your project and compile. These files can be found in the release page.

No external dependencies or additional steps are required.

## Example Usage

```cpp

#include "geomath.hpp"
#include <iostream>

using namespace Arns::Math;

int main() {
    // 1. Geometric Intersections
    Triangle2D triangle({1, 0}, {0, 0}, {0, 1});
    Circle2D circle(Vector2D(-1, 0), 1.0f);

    if (intersect(triangle, circle)) {
        std::cout << "Collision detected!" << std::endl;
    }

    // 2. Method Chaining 
    real_t p = triangle.copy()
                       .scale(2.0f)
                       .perimeter();

    // real_t is the type that is internally used. 
    // Can be set to float or double
    std::cout << "Scaled perimeter: " << p << std::endl;

    // 3. Algebraic Types
    Fraction f1(1, 3);
    Fraction f2(1, 6);
    auto result = f1 + f2; // Returns Fraction(1, 2)
    
    return 0;
}

```

## Status