#include <gtest/gtest.h>
#include <morton.h>
#include <tree.h>

class TestTree : public testing::Test
{
protected:
    void SetUp(){}
    void TearDown(){}
};

TEST_F(TestTree, treeInitialization)
{
    const double BOX[3] = {1.0, 1.0, 1.0};
    ASSERT_NO_FATAL_FAILURE(
            Tree tree = initalizeTree();

            ASSERT_TRUE(tree.leafs != NULL);
            ASSERT_TRUE(tree.firstParticle != NULL);
            ASSERT_TRUE(tree.particleCounts != NULL);

            ASSERT_EQ(0, tree.leafCount);
            ASSERT_EQ(0, tree.particleCounts[0]);

            for (int i = 0; i < MAXLEAFS; ++i) {
                ASSERT_EQ(0u, tree.leafs[i].key);
                ASSERT_EQ(0, tree.firstParticle[i]);
                ASSERT_EQ(0, tree.particleCounts[i]);
            }

            createRootNode(&tree, BOX);

            ASSERT_EQ(1, tree.leafCount);
            ASSERT_EQ(0, tree.particleCounts[0]);

            Morton root = tree.leafs[0];

            ASSERT_EQ(0u, root.key);

            for (int i = 1; i < MAXLEAFS; ++i) {
                ASSERT_EQ(0u, tree.leafs[i].key);
                ASSERT_EQ(0, tree.firstParticle[i]);
                ASSERT_EQ(0, tree.particleCounts[i]);
            }
    );
}

TEST_F(TestTree, particleInRootNode) {
    const double BOX[3] = {1.0, 1.0, 1.0};

    Tree tree = initalizeTree();
    createRootNode(&tree, BOX);

    const double x = drand48(), y = drand48(), z = drand48();

    bool inNode = coordInsideNode(x, y, z, tree.leafs[0], BOX);
    ASSERT_TRUE(inNode);
}
