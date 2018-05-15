#include <gtest/gtest.h>

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    ::testing::GTEST_FLAG (shuffle) = true;

    srand48(time(NULL));

    return RUN_ALL_TESTS();
}



