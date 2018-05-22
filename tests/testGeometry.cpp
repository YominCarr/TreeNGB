#include <gtest/gtest.h>
#include <geometry.h>

class TestGeometry : public testing::Test
{
protected:
    void SetUp(){}
    void TearDown(){}
};

TEST_F(TestGeometry, boxSphereInteraction) {
    double boxCenter[3] = {0.5, 0.5, 0.5};
    double boxSideLength[3] = {0.2, 0.3, 0.25};

    double sphereCener[3] = {0.0, 0.0, 0.0};
    double radius = 0.08;

    ASSERT_FALSE(sphereInsideBox(boxCenter, boxSideLength, sphereCener, radius));

    radius = 0.2;

    ASSERT_FALSE(sphereInsideBox(boxCenter, boxSideLength, sphereCener, radius));

    sphereCener[0] = sphereCener[1] = sphereCener[2] = 0.3;

    ASSERT_TRUE(sphereInsideBox(boxCenter, boxSideLength, sphereCener, radius));

    radius = 0.1;

    ASSERT_TRUE(sphereInsideBox(boxCenter, boxSideLength, sphereCener, radius));

    radius = 0.25;

    ASSERT_FALSE(sphereInsideBox(boxCenter, boxSideLength, sphereCener, radius));

    radius = 0.3;

    ASSERT_FALSE(sphereInsideBox(boxCenter, boxSideLength, sphereCener, radius));

    radius = 0.4;

    ASSERT_FALSE(sphereInsideBox(boxCenter, boxSideLength, sphereCener, radius));
}