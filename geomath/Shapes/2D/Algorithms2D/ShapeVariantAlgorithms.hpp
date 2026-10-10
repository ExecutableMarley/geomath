/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "Algorithms2D.hpp"
#include "../Variants/ShapeVariants.hpp"

namespace Arns::geomath
{

// We can't currently put variant dispatch functions inside Algorithms.hpp

//

//bool isPointInsideShape(const Vector2D& point, const FiniteShapeVariant2D& shell);

//bool isSegmentInsideShape(const Segment2D& segment, const FiniteShapeVariant2D& shell);

bool isShapeInsideShape(const FiniteShapeVariant2D& inner, const FiniteShapeVariant2D& shell);




//

bool intersect(const Ray2D& ray, const FiniteShapeVariant2D& shape, 
    real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D* hitInfo = nullptr);

bool intersect(const Segment2D& segment, const FiniteShapeVariant2D& shape, HitInfo2D* hitInfo = nullptr);

bool intersect(const FiniteShapeVariant2D& shape1, const FiniteShapeVariant2D& shape2);

//

real_t distance(const Vector2D& point, const FiniteShapeVariant2D& shape2, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distance(const Segment2D& segment, const FiniteShapeVariant2D& shape2, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distance(const FiniteShapeVariant2D& shape1, const FiniteShapeVariant2D& shape2, Vector2D* closestPoint1, Vector2D* closestPoint2);

//

SweepResult2D sweep(
    const FiniteShapeVariant2D& shapeA, Vector2D va, 
    const FiniteShapeVariant2D& shapeB, Vector2D vb,
    real_t maxTime = real_t(1), int maxIterations = 15, real_t epsilon = real_t(1e-4));



}

