#include <future>
#include <algorithm>
#include <chrono>
#include <gtest/gtest.h>
#include <morton.h>
#include <tree.h>
#include <particle.h>

#define TEST_TIMEOUT_BEGIN   std::promise<bool> promisedFinished; \
                              auto futureResult = promisedFinished.get_future(); \
                              std::thread([&](std::promise<bool>& finished) {

#define TEST_TIMEOUT_FAIL_END(X)  finished.set_value(true); \
                                   }, std::ref(promisedFinished)).detach(); \
                                   EXPECT_TRUE(futureResult.wait_for(std::chrono::milliseconds(X)) != std::future_status::timeout);

#define TEST_TIMEOUT_SUCCESS_END(X)  finished.set_value(true); \
                                      }, std::ref(promisedFinished)).detach(); \
                                      EXPECT_FALSE(futureResult.wait_for(std::chrono::milliseconds(X)) != std::future_status::timeout);

using namespace std;
using namespace std::chrono;

class TestTree : public testing::Test
{
protected:
    void SetUp(){}
    void TearDown(){}

    unsigned int getRootNode(Tree& tree) {
        unsigned int nodeIndex = 0;
        Morton node = tree.nodes[nodeIndex];
        while (isNotRootNode(node)) {
            nodeIndex = getParentNode(nodeIndex, &tree);
            node = tree.nodes[nodeIndex];
        }
        return nodeIndex;
    }

    int findNGBBruteForce(Particle *P, const int npart, const int ipart, const double hsml, int ngblist[NGBMAX]) {
        int found = 0;
        for (int jpart = 0; jpart < npart; ++jpart) {
            if (ipart == jpart) {
                continue;
            }
            double d = 0.0, d2 = 0.0;
            for (int k = 0; k < 3; ++k) {
                d = P[ipart].Pos[k] - P[jpart].Pos[k];
                d2 += d*d;
            }
            if (d2 <= hsml*hsml) {
                ngblist[found++] = jpart;
            }
        }
        return found;
    }
};

TEST_F(TestTree, treeInitialization)
{
    const double BOX[3] = {1.0, 1.0, 1.0};
    ASSERT_NO_FATAL_FAILURE(
            Tree tree = initalizeTree();

            ASSERT_TRUE(tree.nodes != NULL);
            ASSERT_TRUE(tree.firstParticle != NULL);
            ASSERT_TRUE(tree.particleCounts != NULL);
            ASSERT_TRUE(tree.parentNodes != NULL);
            ASSERT_TRUE(tree.nextNodes != NULL);

            ASSERT_EQ(0, tree.nodeCount);
            ASSERT_EQ(0, tree.particleCounts[0]);

            for (int i = 0; i < MAXLEAVES; ++i) {
                ASSERT_EQ(0u, tree.nodes[i].key);
                ASSERT_EQ(0, tree.firstParticle[i]);
                ASSERT_EQ(0, tree.particleCounts[i]);
                ASSERT_EQ(0u, tree.parentNodes[i]);
                ASSERT_EQ(0u, tree.nextNodes[i]);
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
                ASSERT_EQ(0u, tree.parentNodes[i]);
                ASSERT_EQ(0u, tree.nextNodes[i]);
            }
    );
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
    P[1].leafIndex = 0;
    P[1].Pos[0] = drand48();
    P[1].Pos[1] = drand48();
    P[1].Pos[2] = drand48();

    tree.particleCounts[0] = 1;
    tree.firstParticle[0] = 1;
    for (int i = 1; i < 9; ++i) {
        tree.parentNodes[i] = 1;
    }

    ASSERT_EQ(1, tree.nodeCount);
    ASSERT_EQ(0u, tree.nextNodes[0]);

    ASSERT_NO_FATAL_FAILURE(splitNode(P, 0, &tree, BOX));

    ASSERT_EQ(9, tree.nodeCount);
    ASSERT_EQ(1u, tree.nextNodes[0]);

    const double d[3] = {0.5*BOX[0], 0.5*BOX[1], 0.5*BOX[2]};

    ASSERT_FALSE(nodeIsLeaf(&tree, 0));
    ASSERT_EQ(1u, getFirstSubnodeInNode(0, &tree));

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
                ASSERT_EQ(0u, tree.parentNodes[l]);
                const unsigned int nextNode = (l < 8 ? l+1 : 0);
                ASSERT_EQ(nextNode, tree.nextNodes[l]) << " at l = " << l;
                ASSERT_EQ(nextNode, getNextLeaf(l, &tree)) << " at l = " << l;

                foundParticles += tree.particleCounts[l];
                firstParticles += tree.firstParticle[l];
            }
        }
    }
    ASSERT_EQ(1, foundParticles);
    ASSERT_EQ(1, firstParticles);
}

TEST_F(TestTree, splitNodeTwice) {
    const double BOX[3] = {1.0, 1.0, 1.0};

    Tree tree = initalizeTree();
    createRootNode(&tree, BOX);

    Particle P[2];
    P[1].leafIndex = 0;
    P[1].Pos[0] = drand48();
    P[1].Pos[1] = drand48();
    P[1].Pos[2] = drand48();

    tree.particleCounts[0] = 1;
    tree.firstParticle[0] = 1;
    for (int i = 1; i < 9; ++i) {
        tree.parentNodes[i] = 1;
    }

    ASSERT_EQ(1, tree.nodeCount);
    ASSERT_EQ(0u, tree.nextNodes[0]);

    ASSERT_NO_FATAL_FAILURE(splitNode(P, 0, &tree, BOX));
    ASSERT_NO_FATAL_FAILURE(splitNode(P, 1, &tree, BOX));

    ASSERT_EQ(17, tree.nodeCount);
    ASSERT_EQ(1u, tree.nextNodes[0]);

    const double d[3] = {0.5*BOX[0], 0.5*BOX[1], 0.5*BOX[2]};

    ASSERT_FALSE(nodeIsLeaf(&tree, 0));
    ASSERT_EQ(1u, getFirstSubnodeInNode(0, &tree));

    int foundParticles = 0, firstParticles = 0;

    // First level
    for (int i = 0, l = 1; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k, ++l) {
                double x, y, z;
                key2Coord(tree.nodes[l], &x, &y, &z, BOX);

                ASSERT_EQ(i * d[0], x) << " at l = " << l;
                ASSERT_EQ(j * d[1], y) << " at l = " << l;
                ASSERT_EQ(k * d[2], z) << " at l = " << l;

                ASSERT_EQ(1u, tree.nodes[l].level) << " at l = " << l;

                const unsigned int nextNode = (l == 1 ? 9 : (l < 8 ? l+1 : 0));

                if (l == 1) {
                    ASSERT_FALSE(nodeIsLeaf(&tree, l));
                    ASSERT_EQ(1u, getFirstSubnodeInNode(0, &tree)) << " at l = " << l;
                } else {
                    ASSERT_TRUE(nodeIsLeaf(&tree, l));
                    ASSERT_EQ(nextNode, getNextLeaf(l, &tree)) << " at l = " << l;

                    foundParticles += tree.particleCounts[l];
                    firstParticles += tree.firstParticle[l];
                }

                ASSERT_EQ(0u, tree.parentNodes[l]);
                ASSERT_EQ(nextNode, tree.nextNodes[l]) << " at l = " << l;
            }
        }
    }
    ASSERT_TRUE(foundParticles == 0 || foundParticles == 1);
    ASSERT_TRUE(firstParticles == 0 || firstParticles == 1);

    //Second level
    const double d2[3] = {0.25*BOX[0], 0.25*BOX[1], 0.25*BOX[2]};

    for (int i = 0, l = 9; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k, ++l) {
                double x, y, z;
                key2Coord(tree.nodes[l], &x, &y, &z, BOX);

                ASSERT_EQ(i * d2[0], x) << " at l = " << l;
                ASSERT_EQ(j * d2[1], y) << " at l = " << l;
                ASSERT_EQ(k * d2[2], z) << " at l = " << l;

                ASSERT_EQ(2u, tree.nodes[l].level) << " at l = " << l;

                ASSERT_TRUE(nodeIsLeaf(&tree, l));
                ASSERT_EQ(1u, tree.parentNodes[l]);
                const unsigned int nextNode = (l < 16 ? l+1 : 2); //jump back to leaf in level 1 on the last
                ASSERT_EQ(nextNode, tree.nextNodes[l]) << " at l = " << l;
                ASSERT_EQ(nextNode, getNextLeaf(l, &tree)) << " at l = " << l;

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
    ASSERT_EQ(0u, P[0].leafIndex);
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
    ASSERT_EQ(0u, P[0].leafIndex);
    ASSERT_EQ(0, tree.firstParticle[0]);
    ASSERT_EQ(2, tree.particleCounts[0]);
}

TEST_F(TestTree, rootNodeDetection) {
    const double BOX[3] = {1.0, 1.0, 1.0};
    Tree tree = initalizeTree();
    createRootNode(&tree, BOX);

    ASSERT_FALSE(isNotRootNode(tree.nodes[0]));

    Morton morton;
    morton.key = drand48() * (std::numeric_limits<uint64_t>::max() - 1) + 1;

    ASSERT_TRUE(isNotRootNode(morton));
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

TEST_F(TestTree, leafInsideNode) {
    const double BOX[3] = {1.0, 1.0, 1.0};
    Morton node = coord2Key(0.5, 0.5, 0.5, BOX);
    node.level = 2;

    Morton leaf = coord2Key(0.0, 0.0, 0.0, BOX);
    leaf.level = 3;
    ASSERT_FALSE(nodeContainsLeaf(node, leaf));

    leaf = coord2Key(0.0, 0.6, 0.6, BOX);
    leaf.level = 3;
    ASSERT_FALSE(nodeContainsLeaf(node, leaf));

    leaf = coord2Key(0.5, 0.8, 0.6, BOX);
    leaf.level = 3;
    ASSERT_FALSE(nodeContainsLeaf(node, leaf));

    leaf = coord2Key(0.5, 0.6, 0.6, BOX);
    leaf.level = 3;
    ASSERT_TRUE(nodeContainsLeaf(node, leaf));

    node = coord2Key(0.0, 0.0, 0.0, BOX);
    node.level = 0;

    leaf = coord2Key(0.0, 0.0, 0.25, BOX);
    leaf.level = 2;
    ASSERT_TRUE(nodeContainsLeaf(node, leaf));
}

TEST_F(TestTree, neighbourFindingInLeaf)
{
    const double BOX[3] = {1.0, 1.0, 1.0};

    Particle P[6];
    P[0].Pos[0] = P[0].Pos[1] = P[0].Pos[2] = 0.5; //center of search
    P[1].Pos[0] = P[1].Pos[1] = P[1].Pos[2] = 0.3; //no
    P[2].Pos[0] = P[2].Pos[1] = P[2].Pos[2] = 0.7; //no
    P[3].Pos[0] = P[3].Pos[1] = P[3].Pos[2] = 0.6; //yes
    P[4].Pos[0] = 0.5; P[4].Pos[1] = 0.5; P[4].Pos[2] = 0.7; //yes
    P[5].Pos[0] = 0.4; P[5].Pos[1] = 0.45; P[5].Pos[2] = 0.6; //yes

    // Copy to track back after resort
    Particle Pold[6];
    for (int i = 0; i < 6; ++i) {
        Pold[i] = P[i];
    }

    // Resorts particles by key internally
    Tree tree = buildTree(P, 6, BOX);

    // Track back
    int sort[6];
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            if (Pold[i].Pos[0] == P[j].Pos[0] && Pold[i].Pos[1] == P[j].Pos[1] && Pold[i].Pos[2] == P[j].Pos[2]) {
                sort[i] = j;
                break;
            }
        }
    }

    int ngblist[NGBMAX];
    int found = 0;
    for (int i = 0; i < 3; ++i) {
        ASSERT_NO_FATAL_FAILURE(found += findNeighboursInLeaf(P, sort[0], 0.2, &tree, ngblist, found, P[sort[i]].leafIndex);)
                                    << " at i = " << i;
        ASSERT_EQ(0, found) << " at i = " << i;
    }
    for (int i = 3, j = 0; i < 6; ++i, ++j) {
        ASSERT_NO_FATAL_FAILURE(found += findNeighboursInLeaf(P, sort[0], 0.2, &tree, ngblist, found, P[sort[i]].leafIndex);)
                                    << " i = " << i;
        ASSERT_EQ(j+1, found) << " at i = " << i;
        ASSERT_EQ(sort[i], ngblist[j]) << " at i = " << i;
    }

    freeTreeContents(&tree);
}

TEST_F(TestTree, findRootNode) {
    const double BOX[3] = {1.0, 1.0, 1.0};
    const int N = 100;
    Particle* P = createRandomParticles(N, BOX);

    Tree tree = buildTree(P, N, BOX);

    unsigned int nodeIndex;
    TEST_TIMEOUT_BEGIN
    ASSERT_NO_FATAL_FAILURE(nodeIndex = getRootNode(tree););
    TEST_TIMEOUT_FAIL_END(1000)

    ASSERT_FALSE(isNotRootNode(tree.nodes[nodeIndex]));

    free(P);
}

TEST_F(TestTree, neighbourFindingInNode) {
    const double BOX[3] = {1.0, 1.0, 1.0};

    Particle P[6];
    P[0].Pos[0] = P[0].Pos[1] = P[0].Pos[2] = 0.5; //center of search
    P[1].Pos[0] = P[1].Pos[1] = P[1].Pos[2] = 0.3; //no
    P[2].Pos[0] = P[2].Pos[1] = P[2].Pos[2] = 0.7; //no
    P[3].Pos[0] = P[3].Pos[1] = P[3].Pos[2] = 0.6; //yes
    P[4].Pos[0] = 0.5; P[4].Pos[1] = 0.5; P[4].Pos[2] = 0.7; //yes
    P[5].Pos[0] = 0.4; P[5].Pos[1] = 0.45; P[5].Pos[2] = 0.6; //yes

    // Copy to track back after resort
    Particle Pold[6];
    for (int i = 0; i < 6; ++i) {
        Pold[i] = P[i];
    }

    Tree tree = buildTree(P, 6, BOX);

    // Track back
    int sort[6], isort[6];
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            if (Pold[i].Pos[0] == P[j].Pos[0] && Pold[i].Pos[1] == P[j].Pos[1] && Pold[i].Pos[2] == P[j].Pos[2]) {
                sort[i] = j;
                isort[j] = i;
                break;
            }
        }
    }

    int ngblist[NGBMAX];
    int found;
    unsigned int nodeIndex = getRootNode(tree);
    ASSERT_NO_FATAL_FAILURE(found = findNeighboursInNode(P, sort[0], 0.2, &tree, ngblist, nodeIndex););

    ASSERT_EQ(3, found);

    for (int i = 0; i < 3; ++i) {
        int ngb = isort[ngblist[i]];
        // @todo each only once!
        ASSERT_TRUE(ngb == 3 || ngb == 4 || ngb == 5) << " with ngb = " << ngb;
    }

    freeTreeContents(&tree);
}

TEST_F(TestTree, neighbourFinding) {
    const double BOX[3] = {1.0, 1.0, 1.0};

    Particle P[6];
    P[0].Pos[0] = P[0].Pos[1] = P[0].Pos[2] = 0.5; //center of search
    P[1].Pos[0] = P[1].Pos[1] = P[1].Pos[2] = 0.3; //no
    P[2].Pos[0] = P[2].Pos[1] = P[2].Pos[2] = 0.7; //no
    P[3].Pos[0] = P[3].Pos[1] = P[3].Pos[2] = 0.6; //yes
    P[4].Pos[0] = 0.5; P[4].Pos[1] = 0.5; P[4].Pos[2] = 0.7; //yes
    P[5].Pos[0] = 0.4; P[5].Pos[1] = 0.45; P[5].Pos[2] = 0.6; //yes

    // Copy to track back after resort
    Particle Pold[6];
    for (int i = 0; i < 6; ++i) {
        Pold[i] = P[i];
    }

    Tree tree = buildTree(P, 6, BOX);

    // Track back
    int sort[6], isort[6];
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            if (Pold[i].Pos[0] == P[j].Pos[0] && Pold[i].Pos[1] == P[j].Pos[1] && Pold[i].Pos[2] == P[j].Pos[2]) {
                sort[i] = j;
                isort[j] = i;
                break;
            }
        }
    }

    int ngblist[NGBMAX];
    int found;
    ASSERT_NO_FATAL_FAILURE(found = findNGB(P, sort[0], 0.2, &tree, ngblist, BOX););

    ASSERT_EQ(3, found);

    for (int i = 0; i < 3; ++i) {
        int ngb = isort[ngblist[i]];
        // @todo each only once!
        ASSERT_TRUE(ngb == 3 || ngb == 4 || ngb == 5) << " with ngb = " << ngb;
    }

    freeTreeContents(&tree);
}

TEST_F(TestTree, complexNeighbourFindingAgainstBruteForce) {

    const int NPART = 10000;
    double BOX[3] = {1.0, 1.0, 1.0};

    Particle *P = createRandomParticles(NPART, BOX);

    const int ipart = drand48() * NPART;

    high_resolution_clock::time_point t0 = high_resolution_clock::now();

    Tree tree = buildTree(P, NPART, BOX);

    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    int* ngblistTree = static_cast<int*>(calloc(NGBMAX, sizeof(int)));
    int foundTree = findNGB(P, ipart, 0.1, &tree, ngblistTree, BOX);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();

    int* ngblistBruteForce = static_cast<int*>(calloc(NGBMAX, sizeof(int)));
    int foundBruteForce = findNGBBruteForce(P, NPART, ipart, 0.1, ngblistBruteForce);

    high_resolution_clock::time_point t3 = high_resolution_clock::now();

    ASSERT_EQ(foundBruteForce, foundTree);
    printf("Found %d neighbours with both methods\n", foundTree);

    sort(ngblistTree, ngblistTree+foundTree);
    sort(ngblistBruteForce, ngblistBruteForce+foundBruteForce);

    for (int i = 0; i < foundTree; ++i) {
        ASSERT_EQ(ngblistBruteForce[i], ngblistTree[i]);
    }

    int64_t buildTree = duration_cast<microseconds> ( t1 - t0 ).count();
    int64_t searchTree = duration_cast<microseconds> ( t2 - t1 ).count();
    int64_t searchBruteForce = duration_cast<microseconds> ( t3 - t2 ).count();

    printf("Timings for 1 particle search in a %d particle pool:\
            \n %g ms build tree\n %g ms search tree\n %g ms search brute force\n",
           NPART, buildTree/1000.0, searchTree/1000.0, searchBruteForce/1000.0);

    ASSERT_GE(searchBruteForce, searchTree);
    ASSERT_GE(searchBruteForce*NPART, searchTree*NPART+buildTree);

    freeTreeContents(&tree);
    free(P);
}
