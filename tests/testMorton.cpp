#include <gtest/gtest.h>
#include <morton.h>

class TestMorton : public testing::Test
{
    void SetUp(){}
    void TearDown(){}
};

TEST_F(TestMorton, testCordKeyLoop) {
    const double BOX[3] = {1.0, 1.0, 1.0};
    Morton morton;
    morton.key = drand48() * std::numeric_limits<uint64_t>::max();

    double x, y, z;
    key2Coord(morton, &x, &y, &z, BOX);
    Morton newMorton = coord2Key(x, y, z, BOX);

    EXPECT_EQ(morton.key, newMorton.key);
}