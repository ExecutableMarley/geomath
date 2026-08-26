#include "../../third_party/doctest.h"
#include "Algebra/Matrices/Matrix.hpp"
#include "Algebra/Matrices/Matrix3x3.hpp"
#include "Algebra/Matrices/Matrix4x4.hpp"
#include "Algebra/Matrices/Operators.hpp"

using namespace Arns::Math;

TEST_SUITE("Matrix4x4")
{
    TEST_CASE("Matrix4x4 default constructor")
    {
        Matrix4x4 mat;
        CHECK(mat.rows() == 4);
        CHECK(mat.columns() == 4);

        bool allElementsZero = true;
        for (size_t i = 0; i < mat.rows(); ++i)
            for (size_t j = 0; j < mat.columns(); ++j)
                if (mat(i, j) != real_t{0})
                    allElementsZero = false;
        CHECK(allElementsZero);
    }

    TEST_CASE("Matrix4x4 constructor with initial value")
    {
        Matrix4x4 mat(real_t{5});
        CHECK(mat.rows() == 4);
        CHECK(mat.columns() == 4);

        bool allElementsCorrect = true;
        for (size_t i = 0; i < mat.rows(); ++i)
            for (size_t j = 0; j < mat.columns(); ++j)
                if (mat(i, j) != real_t{5})
                    allElementsCorrect = false;
        CHECK(allElementsCorrect);
    }

    TEST_CASE("Matrix4x4 constructor with individual values")
    {
        Matrix4x4 mat(1, 2, 3, 4,
                      5, 6, 7, 8,
                      9, 10, 11, 12,
                      13, 14, 15, 16);
        CHECK(mat.rows() == 4);
        CHECK(mat.columns() == 4);

        bool allElementsCorrect = true;
        for (size_t i = 0; i < mat.rows(); ++i)
            for (size_t j = 0; j < mat.columns(); ++j)
            {
                const real_t expected = static_cast<real_t>(i * 4 + j + 1);
                if (mat(i, j) != doctest::Approx(expected))
                    allElementsCorrect = false;
            }
        CHECK(allElementsCorrect);
    }

    TEST_CASE("Matrix4x4 copy constructor")
    {
        Matrix4x4 mat1(1, 2, 3, 4,
                       5, 6, 7, 8,
                       9, 10, 11, 12,
                       13, 14, 15, 16);
        Matrix4x4 mat2(mat1);
        CHECK(mat2.rows() == 4);
        CHECK(mat2.columns() == 4);

        bool allElementsCorrect = true;
        for (size_t i = 0; i < mat2.rows(); ++i)
            for (size_t j = 0; j < mat2.columns(); ++j)
                if (mat2(i, j) != doctest::Approx(mat1(i, j)))
                    allElementsCorrect = false;
        CHECK(allElementsCorrect);
    }

    TEST_CASE("Matrix4x4 from IMatrix")
    {
        Matrix mat(4, 4, std::vector<real_t>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16});
        Matrix4x4 mat4x4(mat);
        CHECK(mat4x4.rows() == 4);
        CHECK(mat4x4.columns() == 4);

        bool allElementsCorrect = true;
        for (size_t i = 0; i < mat4x4.rows(); ++i)
            for (size_t j = 0; j < mat4x4.columns(); ++j)
                if (mat4x4(i, j) != doctest::Approx(mat(i, j)))
                    allElementsCorrect = false;
        CHECK(allElementsCorrect);
    }

    TEST_CASE("Matrix4x4 transpose")
    {
        Matrix4x4 mat(1, 2, 3, 4,
                      5, 6, 7, 8,
                      9, 10, 11, 12,
                      13, 14, 15, 16);

        Matrix4x4 expected(1, 5, 9, 13,
                           2, 6, 10, 14,
                           3, 7, 11, 15,
                           4, 8, 12, 16);

        Matrix4x4 transposed = mat.transpose();
        CHECK(transposed == expected);
    }

    TEST_CASE("Matrix4x4 determinant")
    {
        Matrix4x4 mat(1, 2, 3, 4,
                      5, 6, 7, 8,
                      2, 6, 4, 8,
                      3, 1, 1, 2);

        CHECK(mat.determinant() == doctest::Approx(72));
    }

    TEST_CASE("Matrix4x4 Inverse")
    {
        Matrix4x4 mat(3, 4, 5, 6,
                      6, 6, 5, 5,
                      6, 4, 5, 4,
                      3, 4, 5, 3);
        
        REQUIRE(mat.isInvertible());

        Matrix4x4 inv;
        bool invertible = mat.inverse(inv);
        CHECK(invertible);

        // Check that mat * inv is approximately the identity matrix
        bool allElementsCorrect = true;
        Matrix4x4 identity = mat * inv;
        for (size_t i = 0; i < identity.rows(); ++i)
            for (size_t j = 0; j < identity.columns(); ++j)
                if (i == j && identity(i, j) != doctest::Approx(1))
                    allElementsCorrect = false;
                else if (i != j && identity(i, j) != doctest::Approx(0))
                    allElementsCorrect = false;
        CHECK(allElementsCorrect);
    }

    TEST_CASE("Matrix4x4 Inverse of non-invertible matrix")
    {
        Matrix4x4 mat(1, 2, 3, 4,
                      5, 6, 7, 8,
                      9, 10, 11, 12,
                      13, 14, 15, 16); // This matrix is singular (determinant is 0)
        
        CHECK_FALSE(mat.isInvertible());

        Matrix4x4 inv;
        bool invertible = mat.inverse(inv);
        CHECK_FALSE(invertible);
    }

    TEST_CASE("Operator overloads for Matrix4x4")
    {
        SUBCASE("Addition")
        {
            Matrix4x4 mat1(1, 2, 3, 4,
                           5, 6, 7, 8,
                           9, 10, 11, 12,
                           13, 14, 15, 16);
            Matrix4x4 mat2(16, 15, 14, 13,
                           12, 11, 10, 9,
                           8, 7, 6, 5,
                           4, 3, 2, 1);
            Matrix4x4 expected(17, 17, 17, 17,
                               17, 17, 17, 17,
                               17, 17, 17, 17,
                               17, 17, 17, 17);
            CHECK(mat1 + mat2 == expected);
        }

        SUBCASE("Subtraction")
        {
            Matrix4x4 mat1(16, 15, 14, 13,
                           12, 11, 10, 9,
                           8, 7, 6, 5,
                           4, 3, 2, 1);
            Matrix4x4 mat2(1, 2, 3, 4,
                           5, 6, 7, 8,
                           9, 10, 11, 12,
                           13, 14, 15, 16);
            Matrix4x4 expected(15, 13, 11, 9,
                               7, 5, 3, 1,
                               -1, -3, -5, -7,
                               -9, -11, -13, -15);
            CHECK(mat1 - mat2 == expected);
        }

        SUBCASE("Scalar multiplication")
        {
            Matrix4x4 mat(1, 2, 3, 4,
                          5, 6, 7, 8,
                          9, 10, 11, 12,
                          13, 14, 15, 16);
            real_t scalar = 2;
            Matrix4x4 expected(2, 4, 6, 8,
                               10, 12, 14, 16,
                               18, 20, 22, 24,
                               26, 28, 30, 32);
            CHECK(mat * scalar == expected);
        }

        SUBCASE("Matrix multiplication")
        {
            Matrix4x4 mat1(1, 2, 3, 4,
                           5, 6, 7, 8,
                           9, 10, 11, 12,
                           13, 14, 15, 16);
            Matrix4x4 mat2(16, 15, 14, 13,
                           12, 11, 10, 9,
                           8, 7, 6, 5,
                           4, 3, 2, 1);
            Matrix4x4 expected(80, 70, 60, 50,
                               240, 214, 188, 162,
                               400, 358, 316, 274,
                               560, 502, 444, 386);
            CHECK(mat1 * mat2 == expected);
        }

        SUBCASE("Matrix and Vector4D multiplication")
        {
            Matrix4x4 mat(1, 2, 3, 4,
                          5, 6, 7, 8,
                          9, 10, 11, 12,
                          13, 14, 15, 16);
            Vector4D vec(1, 2, 3, 4);
            Vector4D expected(30, 70, 110, 150);
            CHECK(mat * vec == expected);
        }

        SUBCASE("Scalar division")
        {
            Matrix4x4 mat(2, 4, 6, 8,
                          10, 12, 14, 16,
                          18, 20, 22, 24,
                          26, 28, 30, 32);
            real_t scalar = 2;
            Matrix4x4 expected(1, 2, 3, 4,
                               5, 6, 7, 8,
                               9, 10, 11, 12,
                               13, 14, 15, 16);
            CHECK(mat / scalar == expected);
        }
    }

    TEST_CASE("3D homogeneous transforms")
    {
        SUBCASE("Translation")
        {
            Matrix4x4 translationMatrix = Matrix4x4::createTranslation3D(5, 10, 15);
            Vector3D point(2, 3, 4);
            Vector3D transformedPoint = translationMatrix * point;
            CHECK(transformedPoint == Vector3D(7, 13, 19));
        }

        SUBCASE("Scaling")
        {
            Matrix4x4 scalingMatrix = Matrix4x4::createScale3D(2, 3, 4);
            Vector3D point(1, 1, 1);
            Vector3D transformedPoint = scalingMatrix * point;
            CHECK(transformedPoint == Vector3D(2, 3, 4));
        }

        SUBCASE("Rotation around X-axis")
        {
            Matrix4x4 rotationMatrix = Matrix4x4::createRotationXDegs(90);
            Vector3D point(0, 1, 0);
            Vector3D transformedPoint = rotationMatrix * point;
            CHECK(transformedPoint == Vector3D(0, 0, 1));
        }

        SUBCASE("Rotation around Y-axis")
        {
            Matrix4x4 rotationMatrix = Matrix4x4::createRotationYDegs(90);
            Vector3D point(1, 0, 0);
            Vector3D transformedPoint = rotationMatrix * point;
            CHECK(transformedPoint == Vector3D(0, 0, -1));
        }

        SUBCASE("Rotation around Z-axis")
        {
            Matrix4x4 rotationMatrix = Matrix4x4::createRotationZDegs(90);
            Vector3D point(1, 0, 0);
            Vector3D transformedPoint = rotationMatrix * point;
            CHECK(transformedPoint == Vector3D(0, 1, 0));
        }
    }
}