#include <gtest/gtest.h>
#include <particle.h>
#include <morton.h>
#include <tree.h>

class TestParticle : public testing::Test
{
protected:
    void SetUp(){}
    void TearDown(){}
};

TEST_F(TestParticle, particleSwapping) {
    const double BOX[3] = {1.0, 1.0, 1.0};
    const int N = 2;

    Particle* particles = createRandomParticles(N, BOX);
    double p1[3], p2[3];
    for (int i = 0; i < 3; ++i) {
        p1[i] = particles[0].Pos[i];
        p2[i] = particles[1].Pos[i];
    }
    particles[0].leafIndex = 1;
    particles[1].leafIndex = 2;

    swapParticles(particles, 0, 1);

    for (int i = 0; i < 3; ++i) {
        ASSERT_EQ(particles[0].Pos[i], p2[i]);
        ASSERT_EQ(particles[1].Pos[i], p1[i]);
    }
    ASSERT_EQ(particles[0].leafIndex, 2u);
    ASSERT_EQ(particles[1].leafIndex, 1u);

    free(particles);
}

TEST_F(TestParticle, particleSorting) {
    const double BOX[3] = {1.0, 1.0, 1.0};
    const int N = 100;

    Tree tree;
    initalizeTree();

    Particle* particles = createRandomParticles(N, BOX);
    for (int i = 0; i < N; ++i) {
        tree.nodes[i].key = drand48() * std::numeric_limits<uint64_t>::max();
        particles[i].leafIndex = i;
    }

    ASSERT_NO_FATAL_FAILURE(sortParticlesByKey(particles, N, &tree););

    for (int i = 1; i < N; ++i) {
        ASSERT_GE(tree.nodes[particles[i].leafIndex].key, tree.nodes[particles[i-1].leafIndex].key) << " i = " << i;
    }

    freeTreeContents(&tree);
    free(particles);
}