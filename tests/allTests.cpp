#include <gtest/gtest.h>

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    ::testing::GTEST_FLAG (shuffle) = true;

    long seed = time(NULL);
    srand48(seed);
    printf("\nRunning with a random seed of %ld\n\n", seed);

    return RUN_ALL_TESTS();
}
