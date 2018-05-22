#include <gtest/gtest.h>
#include <geometry.h>

class TestGeometry : public testing::Test
{
protected:
    void SetUp(){}
    void TearDown(){}
};

TEST_F(TestGeometry, coordInsideBox)
{
    double boxLowerCoord[3] = {0.3, 0.2, 0.25};
    double boxSideLength[3] = {0.2, 0.3, 0.25};

    double coord[3] = {0.0, 0.0, 0.0};

    ASSERT_FALSE(coordInsideBox(boxLowerCoord, boxSideLength, coord));

    coord[0] = 0.4;

    ASSERT_FALSE(coordInsideBox(boxLowerCoord, boxSideLength, coord));

    coord[1] = 0.5;

    ASSERT_FALSE(coordInsideBox(boxLowerCoord, boxSideLength, coord));

    coord[2] = 0.4;

    ASSERT_TRUE(coordInsideBox(boxLowerCoord, boxSideLength, coord));
}

TEST_F(TestGeometry, fitSphereInBox) {
    double boxLowerCoord[3] = {0.3, 0.2, 0.25};
    double boxSideLength[3] = {0.2, 0.3, 0.25};
    double sphereCenter[3] = {0.4, 0.45, 0.4};

    double maxRadius = getMaxRadiusForSphereInBox(boxLowerCoord, boxSideLength, sphereCenter);

    ASSERT_NEAR(0.05, maxRadius, 0.000001);
}

TEST_F(TestGeometry, sphereInsideBox) {
    double boxLowerCoord[3] = {0.3, 0.2, 0.25};
    double boxSideLength[3] = {0.2, 0.3, 0.25};

    double sphereCenter[3] = {0.0, 0.0, 0.0};
    double radius = 0.2;

    ASSERT_FALSE(sphereInsideBox(boxLowerCoord, boxSideLength, sphereCenter, radius));

    radius = 0.05;

    ASSERT_FALSE(sphereInsideBox(boxLowerCoord, boxSideLength, sphereCenter, radius));

    sphereCenter[0] = sphereCenter[1] = sphereCenter[2] = 0.4;

    ASSERT_TRUE(sphereInsideBox(boxLowerCoord, boxSideLength, sphereCenter, radius));

    radius = 0.0999999;

    ASSERT_TRUE(sphereInsideBox(boxLowerCoord, boxSideLength, sphereCenter, radius));

    radius = 0.25;

    ASSERT_FALSE(sphereInsideBox(boxLowerCoord, boxSideLength, sphereCenter, radius));

    radius = 0.3;

    ASSERT_FALSE(sphereInsideBox(boxLowerCoord, boxSideLength, sphereCenter, radius));

    radius = 0.4;

    ASSERT_FALSE(sphereInsideBox(boxLowerCoord, boxSideLength, sphereCenter, radius));
}