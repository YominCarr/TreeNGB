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

            ASSERT_TRUE(tree.nodes != NULL);
            ASSERT_TRUE(tree.firstParticle != NULL);
            ASSERT_TRUE(tree.particleCounts != NULL);

            ASSERT_EQ(0, tree.nodeCount);
            ASSERT_EQ(0, tree.particleCounts[0]);

            for (int i = 0; i < MAXLEAVES; ++i) {
                ASSERT_EQ(0u, tree.nodes[i].key);
                ASSERT_EQ(0, tree.firstParticle[i]);
                ASSERT_EQ(0, tree.particleCounts[i]);
            }

            createRootNode(&tree, BOX);

            ASSERT_EQ(1, tree.nodeCount);
            ASSERT_EQ(0, tree.particleCounts[0]);

            Morton root = tree.nodes[0];

            ASSERT_EQ(0u, root.key);

            for (int i = 1; i < MAXLEAVES; ++i) {
                ASSERT_EQ(0u, tree.nodes[i].key);
                ASSERT_EQ(0, tree.firstParticle[i]);
                ASSERT_EQ(0, tree.particleCounts[i]);
            }
    );
}

TEST_F(TestTree, nodeSizeCalculation) {
    const double BOX[3] = {1.0, 1.0, 1.0};
    double sideLength[3];
    Morton key;

    key.level = 0;
    getNodeSize(sideLength, key, BOX);
    for (int i = 0; i < 3; ++i) {
        ASSERT_EQ(BOX[i], sideLength[i]) << " at i = " << i;
    }

    key.level = 1;
    getNodeSize(sideLength, key, BOX);
    for (int i = 0; i < 3; ++i) {
        ASSERT_EQ(0.5*BOX[i], sideLength[i]) << " at i = " << i;
    }
}

TEST_F(TestTree, particleInRootNode) {
    const double BOX[3] = {1.0, 1.0, 1.0};

    Tree tree = initalizeTree();
    createRootNode(&tree, BOX);

    const double x = drand48(), y = drand48(), z = drand48();

    bool inNode = coordInsideNode(x, y, z, tree.nodes[0], BOX);
    ASSERT_TRUE(inNode);
}

TEST_F(TestTree, particleInNode) {
    const double BOX[3] = {1.0, 1.0, 1.0};

    Morton node = coord2Key(0.0, 0.5, 0.0, BOX);
    node.level = 1;

    const double x = 0.2, y = 0.7, z = 0.3;

    bool inNode = coordInsideNode(x, y, z, node, BOX);
    ASSERT_TRUE(inNode);
}

TEST_F(TestTree, splitRootNode) {
    const double BOX[3] = {1.0, 1.0, 1.0};

    Tree tree = initalizeTree();
    createRootNode(&tree, BOX);

    Particle P[2];
    P[1].leaf = tree.nodes[0];
    P[1].Pos[0] = drand48();
    P[1].Pos[1] = drand48();
    P[1].Pos[2] = drand48();

    tree.particleCounts[0] = 1;
    tree.firstParticle[0] = 1;

    ASSERT_EQ(1, tree.nodeCount);

    ASSERT_NO_FATAL_FAILURE(splitNode(P, 0, &tree, BOX));

    ASSERT_EQ(9, tree.nodeCount);

    const double d[3] = {0.5*BOX[0], 0.5*BOX[1], 0.5*BOX[2]};

    ASSERT_FALSE(nodeIsLeaf(&tree, 0));

    int foundParticles = 0, firstParticles = 0;
    for (int i = 0, l = 1; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k, ++l) {
                double x, y, z;
                key2Coord(tree.nodes[l], &x, &y, &z, BOX);

                ASSERT_EQ(i * d[0], x) << " at l = " << l;
                ASSERT_EQ(j * d[1], y) << " at l = " << l;
                ASSERT_EQ(k * d[2], z) << " at l = " << l;

                ASSERT_EQ(1u, tree.nodes[l].level) << " at l = " << l;

                ASSERT_TRUE(nodeIsLeaf(&tree, l));

                foundParticles += tree.particleCounts[l];
                firstParticles += tree.firstParticle[l];
            }
        }
    }
    ASSERT_EQ(1, foundParticles);
    ASSERT_EQ(1, firstParticles);
}

TEST_F(TestTree, assignParticleToTree) {
    const double BOX[3] = {1.0, 1.0, 1.0};

    Tree tree = initalizeTree();
    createRootNode(&tree, BOX);
    tree.firstParticle[0] = INT32_MAX;

    Particle P[1];
    P[0].Pos[0] = drand48();
    P[0].Pos[1] = drand48();
    P[0].Pos[2] = drand48();

    unsigned int leaf = UINT32_MAX;
    ASSERT_NO_FATAL_FAILURE(leaf = assignParticleToTree(P, 0, &tree, BOX));

    ASSERT_EQ(0u, leaf);
    ASSERT_EQ(tree.nodes[0].key, P[0].leaf.key);
    ASSERT_EQ(0, tree.firstParticle[0]);
    ASSERT_EQ(1, tree.particleCounts[0]);
}

TEST_F(TestTree, assignParticleToOccupiedNode) {
    const double BOX[3] = {1.0, 1.0, 1.0};

    Tree tree = initalizeTree();
    createRootNode(&tree, BOX);
    tree.firstParticle[0] = INT32_MAX;

    Particle P[2];
    for (int i = 0; i < 2; ++i) {
        P[i].Pos[0] = drand48();
        P[i].Pos[1] = drand48();
        P[i].Pos[2] = drand48();
    }

    int leaf = -1;
    ASSERT_NO_FATAL_FAILURE(assignParticleToTree(P, 0, &tree, BOX));
    ASSERT_NO_FATAL_FAILURE(leaf = assignParticleToTree(P, 1, &tree, BOX));

    ASSERT_EQ(0, leaf);
    ASSERT_EQ(tree.nodes[0].key, P[0].leaf.key);
    ASSERT_EQ(0, tree.firstParticle[0]);
    ASSERT_EQ(2, tree.particleCounts[0]);
}

TEST_F(TestTree, rootNodeDetection) {
    const double BOX[3] = {1.0, 1.0, 1.0};
    Tree tree = initalizeTree();
    createRootNode(&tree, BOX);

    ASSERT_TRUE(isNotRootNode(tree.nodes[0]));

    Morton morton;
    morton.key = drand48() * std::numeric_limits<uint64_t>::max();

    ASSERT_FALSE(isNotRootNode(morton));
}

TEST_F(TestTree, nodeToBoxTransformation) {
    const double BOX[3] = {1.0, 1.0, 1.0};
    double lowerCoords[3], sideLength[3];

    Morton node = coord2Key(0.25, 0.5, 0.75, BOX);
    node.level = 2;

    ASSERT_NO_FATAL_FAILURE(nodeToBox(node, lowerCoords, sideLength, BOX));

    for (int i = 0; i < 3; ++i) {
        ASSERT_EQ(0.25, sideLength[i]) << " at i = " << i;
    }
    ASSERT_EQ(0.25, lowerCoords[0]);
    ASSERT_EQ(0.5, lowerCoords[1]);
    ASSERT_EQ(0.75, lowerCoords[2]);
}

TEST_F(TestTree, nodeSphereInteraction) {
    const double BOX[3] = {1.0, 1.0, 1.0};
    Morton node = coord2Key(0.25, 0.5, 0.75, BOX);
    node.level = 2;
    const double center[3] = {0.4, 0.65, 0.9};
    double radius = 0.05;

    ASSERT_TRUE(nodeSurroundsSphere(node, center, radius, BOX));

    radius = 0.2;

    ASSERT_FALSE(nodeSurroundsSphere(node, center, radius, BOX));
}
