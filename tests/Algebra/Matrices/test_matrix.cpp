#include "../../third_party/doctest.h"
#include "Algebra/Matrices/Matrix.hpp"
#include "Algebra/Matrices/Matrix3x3.hpp"
#include "Algebra/Matrices/Matrix4x4.hpp"
#include "Algebra/Matrices/Operators.hpp"

using namespace Arns::Math;

TEST_SUITE("Matrix")
{
    TEST_CASE("Matrix default constructor")
    {
        Matrix mat;
        CHECK(mat.rows() == 0);
        CHECK(mat.columns() == 0);
    }

    TEST_CASE("Matrix constructor with dimensions and initial value")
    {
        Matrix mat(2, 3, real_t{1});
        CHECK(mat.rows() == 2);
        CHECK(mat.columns() == 3);

        bool allElementsCorrect = true;
        for (size_t i = 0; i < mat.rows(); ++i)
            for (size_t j = 0; j < mat.columns(); ++j)
                if (mat(i, j) != real_t{1})
                    allElementsCorrect = false;
        CHECK(allElementsCorrect);
    }

    TEST_CASE("Matrix constructor with dimensions and data vector")
    {
        std::vector<real_t> data = {1, 2, 3, 4, 5, 6};
        Matrix mat(2, 3, data);
        CHECK(mat.rows() == 2);
        CHECK(mat.columns() == 3);

        bool allElementsCorrect = true;
        for (size_t i = 0; i < mat.rows(); ++i)
            for (size_t j = 0; j < mat.columns(); ++j)
                if (mat(i, j) != data[i * mat.columns() + j])
                    allElementsCorrect = false;
        CHECK(allElementsCorrect);
    }

    TEST_CASE("Matrix copy constructor")
    {
        Matrix mat1(2, 2, 5.0f);
        Matrix mat2(mat1);
        CHECK(mat2.rows() == 2);
        CHECK(mat2.columns() == 2);
        for (size_t i = 0; i < mat2.rows(); ++i)
            for (size_t j = 0; j < mat2.columns(); ++j)
                CHECK(mat2(i, j) == doctest::Approx(5.0f));
    }

    TEST_CASE("Matrix from IMatrix")
    {
        Matrix3x3 m3x3(1, 2, 3,
                       4, 5, 6,
                       7, 8, 9);
        Matrix m(m3x3);
        CHECK(m.rows() == 3);
        CHECK(m.columns() == 3);
        for (size_t i = 0; i < m.rows(); ++i)
            for (size_t j = 0; j < m.columns(); ++j)
                CHECK(m(i, j) == doctest::Approx(m3x3(i, j)));
    }

    TEST_CASE("Matrix transpose")
    {
        Matrix mat(2, 3);
        mat(0, 0) = 1; mat(0, 1) = 2; mat(0, 2) = 3;
        mat(1, 0) = 4; mat(1, 1) = 5; mat(1, 2) = 6;

        Matrix expected(3, 2);
        expected(0, 0) = 1; expected(0, 1) = 4;
        expected(1, 0) = 2; expected(1, 1) = 5;
        expected(2, 0) = 3; expected(2, 1) = 6;

        Matrix transposed = mat.transpose();
        CHECK(transposed == expected);
    }

    TEST_CASE("Operator overloads")
    {
        SUBCASE("Equality")
        {
            Matrix mat1{{1, 2, 3},
                        {4, 5, 6},
                        {7, 8, 9}};
            Matrix mat2{{1, 2, 3},
                        {4, 5, 6},
                        {7, 8, 9}};
            CHECK(mat1 == mat2);
        }

        SUBCASE("Inequality")
        {
            Matrix mat1{{1, 2, 3},
                        {4, 5, 6},
                        {7, 8, 9}};
            Matrix mat2{{9, 8, 7},
                        {6, 5, 4},
                        {3, 2, 1}};
            CHECK(mat1 != mat2);
        }

        SUBCASE("Addition")
        {
            Matrix mat1 = {{1, 2, 3},
                           {4, 5, 6},
                           {7, 8, 9}};
            Matrix mat2 = {{9, 8, 7},
                           {6, 5, 4},
                           {3, 2, 1}};

            Matrix expected = {{10, 10, 10},
                               {10, 10, 10},
                               {10, 10, 10}};
            CHECK(mat1 + mat2 == expected);
        }

        SUBCASE("Scalar addition")
        {
            Matrix mat = {{1, 2, 3},
                           {4, 5, 6},
                           {7, 8, 9}};
            real_t scalar = 1;

            Matrix expected = {{2, 3, 4},
                               {5, 6, 7},
                               {8, 9, 10}};
            CHECK(mat + scalar == expected);
        }

        SUBCASE("Subtraction")
        {
            Matrix mat1 = {{9, 8, 7},
                           {6, 5, 4},
                           {3, 2, 1}};
            Matrix mat2 = {{1, 2, 3},
                           {4, 5, 6},
                           {7, 8, 9}};

            Matrix expected = {{8, 6, 4},
                               {2, 0, -2},
                               {-4, -6, -8}};
            CHECK(mat1 - mat2 == expected);
        }

        SUBCASE("Scalar subtraction")
        {
            Matrix mat = {{9, 8, 7},
                           {6, 5, 4},
                           {3, 2, 1}};
            real_t scalar = 1;

            Matrix expected = {{8, 7, 6},
                               {5, 4, 3},
                               {2, 1, 0}};
            CHECK(mat - scalar == expected);
        }

        SUBCASE("Matrix multiplication")
        {
            Matrix mat1 = {{1, 2, 3},
                            {4, 5, 6},
                            {7, 8, 9}};
            Matrix mat2 = {{9, 8, 7},
                            {6, 5, 4},
                            {3, 2, 1}};

            Matrix expected = {{30, 24, 18},
                               {84, 69, 54},
                               {138, 114, 90}};
            CHECK(mat1 * mat2 == expected);
        }

        SUBCASE("Scalar multiplication")
        {
            Matrix mat = {{1, 2, 3},
                           {4, 5, 6},
                           {7, 8, 9}};
            real_t scalar = 2;

            Matrix expected = {{2, 4, 6},
                               {8, 10, 12},
                               {14, 16, 18}};
            CHECK(mat * scalar == expected);
        }

        /* Todo: 
        SUBCASE("Matrix and Vector multiplication")
        {
            Matrix mat = {{1, 2, 3},
                           {4, 5, 6},
                           {7, 8, 9}};
            Vector3D vec(1, 2, 3);

            Vector3D expected(14, 32, 50);
            CHECK(mat * vec == expected);
        }*/

        SUBCASE("Scalar division")
        {
            Matrix mat = {{2, 4, 6},
                           {8, 10, 12},
                           {14, 16, 18}};
            real_t scalar = 2;

            Matrix expected = {{1, 2, 3},
                               {4, 5, 6},
                               {7, 8, 9}};
            CHECK(mat / scalar == expected);
        }
    }
}