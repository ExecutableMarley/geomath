#include "../../third_party/doctest.h"
#include "Algebra/Matrices/Matrix.hpp"
#include "Algebra/Matrices/Matrix3x3.hpp"
#include "Algebra/Matrices/Matrix4x4.hpp"
#include "Algebra/Matrices/Operators.hpp"

using namespace Arns::Math;

TEST_SUITE("Matrix3x3")
{
    TEST_CASE("Matrix3x3 default constructor")
    {
        Matrix3x3 mat;
        CHECK(mat.rows() == 3);
        CHECK(mat.columns() == 3);

        bool allElementsZero = true;
        for (size_t i = 0; i < mat.rows(); ++i)
            for (size_t j = 0; j < mat.columns(); ++j)
                if (mat(i, j) != real_t{0})
                    allElementsZero = false;
        CHECK(allElementsZero);
    }

    TEST_CASE("Matrix3x3 constructor with initial value")
    {
        Matrix3x3 mat(real_t{5});
        CHECK(mat.rows() == 3);
        CHECK(mat.columns() == 3);

        bool allElementsCorrect = true;
        for (size_t i = 0; i < mat.rows(); ++i)
            for (size_t j = 0; j < mat.columns(); ++j)
                if (mat(i, j) != real_t{5})
                    allElementsCorrect = false;
        CHECK(allElementsCorrect);
    }

    TEST_CASE("Matrix3x3 constructor with individual values")
    {
        Matrix3x3 mat(1, 2, 3,
                      4, 5, 6,
                      7, 8, 9);
        CHECK(mat.rows() == 3);
        CHECK(mat.columns() == 3);

        bool allElementsCorrect = true;
        for (size_t i = 0; i < mat.rows(); ++i)
            for (size_t j = 0; j < mat.columns(); ++j)
            {
                if (mat(i, j) != (real_t)i * (real_t)3 + (real_t)j + (real_t)1)
                    allElementsCorrect = false;
            }  
        CHECK(allElementsCorrect);
    }

    TEST_CASE("Matrix3x3 copy constructor")
    {
        Matrix3x3 mat1(1, 2, 3,
                       4, 5, 6,
                       7, 8, 9);
        Matrix3x3 mat2(mat1);
        CHECK(mat2.rows() == 3);
        CHECK(mat2.columns() == 3);

        bool allElementsCorrect = true;
        for (size_t i = 0; i < mat2.rows(); ++i)
            for (size_t j = 0; j < mat2.columns(); ++j)
                if (mat2(i, j) != doctest::Approx(mat1(i, j)))
                    allElementsCorrect = false;
        CHECK(allElementsCorrect);
    }

    TEST_CASE("Matrix3x3 from IMatrix")
    {
        Matrix mat(3, 3, std::vector<real_t>{1, 2, 3, 4, 5, 6, 7, 8, 9});
        Matrix3x3 mat3x3(mat);
        CHECK(mat3x3.rows() == 3);
        CHECK(mat3x3.columns() == 3);

        bool allElementsCorrect = true;
        for (size_t i = 0; i < mat3x3.rows(); ++i)
            for (size_t j = 0; j < mat3x3.columns(); ++j)
                if (mat3x3(i, j) != doctest::Approx(mat(i, j)))
                    allElementsCorrect = false;
        CHECK(allElementsCorrect);
    }

    TEST_CASE("Matrix3x3 transpose")
    {
        Matrix3x3 mat(1, 2, 3,
                      4, 5, 6,
                      7, 8, 9);
        
        Matrix3x3 expected(1, 4, 7,
                          2, 5, 8,
                          3, 6, 9);
        
        Matrix3x3 transposed = mat.transpose();
        CHECK(mat.transpose() == expected);
    }

    TEST_CASE("Matrix3x3 determinant")
    {
        Matrix3x3 mat(1, 2, 3,
                      0, 1, 4,
                      5, 6, 0);
        
        real_t det = mat.determinant();
        CHECK(det == doctest::Approx(1 * (1 * 0 - 4 * 6) - 
                                      2 * (0 * 0 - 4 * 5) + 
                                      3 * (0 * 6 - 1 * 5)));
    }

    TEST_CASE("Matrix3x3 Inverse")
    {
        Matrix3x3 mat(1, 2, 3,
                      0, 1, 4,
                      5, 6, 0);
        
        REQUIRE(mat.isInvertible());

        Matrix3x3 inv;
        bool invertible = mat.inverse(inv);
        CHECK(invertible);

        // Check that mat * inv is approximately the identity matrix
        CHECK(mat * inv == Matrix3x3::identity());
    }

    TEST_CASE("Matrix3x3 Inverse of non-invertible matrix")
    {
        Matrix3x3 mat(1, 2, 3,
                      4, 5, 6,
                      7, 8, 9); // This matrix is singular (determinant is 0)
        
        CHECK_FALSE(mat.isInvertible());

        Matrix3x3 inv;
        bool invertible = mat.inverse(inv);
        CHECK_FALSE(invertible);
    }

    TEST_CASE("Operator overloads for Matrix3x3")
    {
        SUBCASE("Addition")
        {
            Matrix3x3 mat1(1, 2, 3,
                           4, 5, 6,
                           7, 8, 9);
            Matrix3x3 mat2(9, 8, 7,
                           6, 5, 4,
                           3, 2, 1);
            Matrix3x3 expected(10, 10, 10,
                               10, 10, 10,
                               10, 10, 10);
            CHECK(mat1 + mat2 == expected);
        }

        SUBCASE("Subtraction")
        {
            Matrix3x3 mat1(9, 8, 7,
                           6, 5, 4,
                           3, 2, 1);
            Matrix3x3 mat2(1, 2, 3,
                           4, 5, 6,
                           7, 8, 9);
            Matrix3x3 expected(8, 6, 4,
                               2, 0, -2,
                               -4, -6, -8);
            CHECK(mat1 - mat2 == expected);
        }

        SUBCASE("Scalar multiplication")
        {
            Matrix3x3 mat(1, 2, 3,
                          4, 5, 6,
                          7, 8, 9);
            real_t scalar = 2;
            Matrix3x3 expected(2, 4, 6,
                               8, 10, 12,
                               14, 16, 18);
            CHECK(mat * scalar == expected);
        }

        SUBCASE("Matrix multiplication")
        {
            Matrix3x3 mat1(1, 2, 3,
                           4, 5, 6,
                           7, 8, 9);
            Matrix3x3 mat2(9, 8, 7,
                           6, 5, 4,
                           3, 2, 1);
            Matrix3x3 expected(30, 24, 18,
                               84, 69, 54,
                               138, 114, 90);
            CHECK(mat1 * mat2 == expected);
        }

        SUBCASE("Matrix and Vector multiplication")
        {
            Matrix3x3 mat(1, 2, 3,
                          4, 5, 6,
                          7, 8, 9);
            Vector3D vec(1, 2, 3);
            Vector3D expected(14, 32, 50);
            CHECK(mat * vec == expected);
        }

        SUBCASE("Scalar division")
        {
            Matrix3x3 mat(2, 4, 6,
                          8, 10, 12,
                          14, 16, 18);
            real_t scalar = 2;
            Matrix3x3 expected(1, 2, 3,
                               4, 5, 6,
                               7, 8, 9);
            CHECK(mat / scalar == expected);
        }
    }

    TEST_CASE("2D homogeneous transforms")
    {
        SUBCASE("Translation")
        {
            Matrix3x3 translationMatrix = Matrix3x3::createTranslation2D(Vector2D(5, 10));
            Vector2D point(2, 3);
            Vector2D transformedPoint = translationMatrix * point;
            CHECK(transformedPoint == Vector2D(7, 13));
        }

        SUBCASE("Scaling")
        {
            Matrix3x3 scalingMatrix = Matrix3x3::createScale2D(2, 3);
            Vector2D point(4, 5);
            Vector2D transformedPoint = scalingMatrix * point;
            CHECK(transformedPoint == Vector2D(8, 15));
        }

        SUBCASE("Rotation")
        {
            Matrix3x3 rotationMatrix = Matrix3x3::createRotationDegs2D(90); // 90 degrees
            Vector2D point(1, 0);
            Vector2D transformedPoint = rotationMatrix * point;
            CHECK(transformedPoint == Vector2D(0, 1));
        }

        SUBCASE("Combined Transform")
        {
            Matrix3x3 translationMatrix = Matrix3x3::createTranslation2D(Vector2D(5, 10));
            Matrix3x3 scalingMatrix = Matrix3x3::createScale2D(2, 3);
            Matrix3x3 rotationMatrix = Matrix3x3::createRotationDegs2D(90); // 90 degrees

            // Combined transformation: first scale, then rotate, then translate
            Matrix3x3 combinedTransform = translationMatrix * rotationMatrix * scalingMatrix;

            Vector2D point(1, 1);
            Vector2D transformedPoint = combinedTransform * point;

            // Expected result after scaling (2, 3), rotating 90 degrees, and translating (5, 10)
            CHECK(transformedPoint == Vector2D(2, 12));
        }
    }

    TEST_CASE("3D Linear transforms")
    {
        SUBCASE("Scaling")
        {
            Matrix3x3 scalingMatrix = Matrix3x3::createScale3D(2, 3, 4);
            Vector3D point(1, 1, 1);
            Vector3D transformedPoint = scalingMatrix * point;
            CHECK(transformedPoint == Vector3D(2, 3, 4));
        }

        SUBCASE("Rotation around X-axis")
        {
            Matrix3x3 rotationMatrix = Matrix3x3::createRotationXDegs(90); // 90 degrees around X-axis
            Vector3D point(0, 1, 0);
            Vector3D transformedPoint = rotationMatrix * point;
            CHECK(transformedPoint == Vector3D(0, 0, 1));
        }

        SUBCASE("Rotation around Y-axis")
        {
            Matrix3x3 rotationMatrix = Matrix3x3::createRotationYDegs(90); // 90 degrees around Y-axis
            Vector3D point(1, 0, 0);
            Vector3D transformedPoint = rotationMatrix * point;
            CHECK(transformedPoint == Vector3D(0, 0, -1));
        }

        SUBCASE("Rotation around Z-axis")
        {
            Matrix3x3 rotationMatrix = Matrix3x3::createRotationZDegs(90); // 90 degrees around Z-axis
            Vector3D point(1, 0, 0);
            Vector3D transformedPoint = rotationMatrix * point;
            CHECK(transformedPoint == Vector3D(0, 1, 0));
        }
    }
}