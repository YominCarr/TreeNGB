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
    EXPECT_EQ(0, morton.x) << stage;
    EXPECT_EQ(0, morton.y) << stage;
    EXPECT_EQ(0, morton.z) << stage;
    EXPECT_EQ(0, morton.level) << stage;
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

TEST_F(TestMorton, structAssignement) {
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

TEST_F(TestMorton, cordKeyLoop) {
    const double BOX[3] = {1.0, 1.0, 1.0};
    Morton morton;
    morton.key = drand48() * std::numeric_limits<uint64_t>::max();

    double x, y, z;
    key2Coord(morton, &x, &y, &z, BOX);
    Morton newMorton = coord2Key(x, y, z, BOX);

    ASSERT_EQ(morton.key, newMorton.key);
}