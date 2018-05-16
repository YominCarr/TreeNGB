#include <gtest/gtest.h>
#include <morton.h>

class TestMorton : public testing::Test
{
protected:
    void SetUp(){}
    void TearDown(){}

    void expectZero(Morton& morton, const char* stage);
};

void TestMorton::expectZero(Morton& morton, const char* stage) {
    EXPECT_EQ(0u, morton.key) << stage;
    EXPECT_EQ(0u, morton.x) << stage;
    EXPECT_EQ(0u, morton.y) << stage;
    EXPECT_EQ(0u, morton.z) << stage;
    EXPECT_EQ(0u, morton.level) << stage;
}

TEST_F(TestMorton, keySize) {
    EXPECT_EQ(64u, 8*sizeof(Morton));

    Morton morton;

    EXPECT_EQ(64u, 8*sizeof(morton.key));
}

TEST_F(TestMorton, translateCoordDoubleLoop) {
    const double box = 1.0;

    double d = drand48();
    COORD c = translateCoordFromDouble(d, box);
    double d2 = translateCoordToDouble(c, box);

    ASSERT_NEAR(d, d2, 0.0001*d);
}

TEST_F(TestMorton, zeroStructAssignement) {
    Morton morton;
    //Not zeroed after construnction!
    //this->expectZero(morton, "after constructor");

    morton.key = drand48() * std::numeric_limits<uint64_t>::max();
    morton.key = 0;
    this->expectZero(morton, "after key assignment");

    morton.key = drand48() * std::numeric_limits<uint64_t>::max();
    morton.x = 0;
    morton.y = 0;
    morton.z = 0;
    morton.level = 0;
    this->expectZero(morton, "after bitfield assignment");
}

TEST_F(TestMorton, halfSizeNode2Key) {
    const double BOX[3] = {1.0, 1.0, 1.0};
    const double x = 0.0, y = 0.0, z = 0.5 * BOX[2];

    Morton node = coord2Key(x, y, z, BOX);
    double x2, y2, z2;
    key2Coord(node, &x2, &y2, &z2, BOX);

    ASSERT_EQ(x, x2);
    ASSERT_EQ(y, y2);
    ASSERT_EQ(z, z2);

    const double level = 1;
    node.level = level;

    ASSERT_EQ(x, x2);
    ASSERT_EQ(y, y2);
    ASSERT_EQ(z, z2);
    ASSERT_EQ(level, node.level);
}

TEST_F(TestMorton, keyCoordLoop) {
    const double BOX[3] = {1.0, 1.0, 1.0};
    Morton morton;
    morton.key = drand48() * std::numeric_limits<uint64_t>::max();

    double x, y, z;
    key2Coord(morton, &x, &y, &z, BOX);
    Morton newMorton = coord2Key(x, y, z, BOX);

    ASSERT_EQ(morton.key, newMorton.key);
}

TEST_F(TestMorton, coordKeyLoop) {
    const double BOX[3] = {1.0, 1.0, 1.0};
    const double x = drand48() * BOX[0];
    const double y = drand48() * BOX[1];
    const double z = drand48() * BOX[2];

    Morton morton = coord2Key(x, y, z, BOX);

    double x2, y2, z2;
    key2Coord(morton, &x2, &y2, &z2, BOX);

    ASSERT_NEAR(x, x2, 0.0001 * x);
    ASSERT_NEAR(y, y2, 0.0001 * y);
    ASSERT_NEAR(z, z2, 0.0001 * z);
}