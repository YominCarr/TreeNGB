#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "tree.h"
#include "geometry.h"
#include "morton.h"

Tree buildTree(Particle *P, const int npart, const double BOX[3]) {
    fprintf(stderr, "Implement a parallel tree build\n");

    Tree tree = initalizeTree();

    createRootNode(&tree, BOX);
    buildTreeSerial(P, npart, &tree, BOX);

    return tree;
}

Tree initalizeTree() {
    Tree tree;

    tree.nodes = calloc(MAXLEAVES, sizeof(Morton));
    tree.firstParticle = calloc(MAXLEAVES, sizeof(int));
    tree.particleCounts = calloc(MAXLEAVES, sizeof(int));
    tree.parentNodes = calloc(MAXLEAVES, sizeof(unsigned int));
    tree.nextNodes = calloc(MAXLEAVES, sizeof(unsigned int));

    tree.nodeCount = 0;

    return tree;
}

void createRootNode(Tree *tree, const double *BOX) {
    tree->nodeCount = 1;

    tree->nodes[0] = coord2Key(0.0, 0.0, 0.0, BOX);
    tree->nodes[0].level = 0;

}

void buildTreeSerial(Particle *P, const int npart, Tree *tree, const double *BOX) {
    for (int ipart = 0; ipart < npart; ++ipart) {
        unsigned int leaf = assignParticleToTree(P, ipart, tree, BOX);

        while (tree->particleCounts[leaf] > MAXLEAFSIZE) {
            if (! treeHasSpaceForSplittingOnce(tree, tree->nodes[leaf].level+1)) {
                fprintf(stderr, "Trying to refine tree beyond capacity (leaf count = %d; new depth = %d), aborting...",
                        tree->nodeCount, tree->nodes[leaf].level+1);
                exit(1);
            }

            splitNode(P, leaf, tree, BOX);
            leaf = assignParticleToTree(P, ipart, tree, BOX);
        }
    }

    sortParticlesByKey(P, npart, tree);
    //After resorting particles we need to readjust the tree
    // @todo maybe absorb this directly into the resorting?
    setParticleRangesInTree(P, npart, tree);
}

bool treeHasSpaceForSplittingOnce(Tree* tree, int newDepth)
{
    return tree->nodeCount + 8 <= MAXLEAVES && newDepth <= MAXDEPTH;
}

void splitNode(Particle *P, const unsigned int l, Tree *tree, const double BOX[3]) {
    const Morton parent = tree->nodes[l];
    double pX, pY, pZ; //Assume these are the corner with the smallest coord
    key2Coord(parent, &pX, &pY, &pZ, BOX);
    const int parentLevel = key2Depth(parent);

    //Particle to redistribute
    int epart = -1;
    if (tree->particleCounts[l] > 0) {
        epart = tree->firstParticle[l];
    }

    // @todo Do I need to track these stats for the non leaf nodes?
    tree->firstParticle[l] = -1; // Minus to mark it is not a leaf
    tree->particleCounts[l] = 0;

    const unsigned int nextLeaf = tree->nextNodes[l];
    tree->nextNodes[l] = tree->nodeCount; //First subnode

    double newSize[3];
    for (int i = 0; i < 3; ++i) {
        newSize[i] = BOX[i] / (2 << parentLevel);
    }

    for (int i = 0, s = tree->nodeCount; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k, ++s) {
                //! @todo can be improved by working solely with the integers
                const double x = pX + i * newSize[0];
                const double y = pY + j * newSize[1];
                const double z = pZ + k * newSize[2];
                Morton node = coord2Key(x, y, z, BOX);
                node.level = parentLevel + 1;

                tree->nodes[s] = node;
                tree->parentNodes[s] = l;
                // @todo do that here or later?
                if (s > tree->nodeCount) {
                    tree->nextNodes[s-1] = s;
                }

                if (epart >= 0 && coordInsideNode(P[epart].Pos[0], P[epart].Pos[1], P[epart].Pos[2], node, BOX)) {
                    P[epart].leafIndex = s;
                    tree->particleCounts[s] = 1;
                    tree->firstParticle[s] = epart;
                }
            }
        }
    }
    tree->nextNodes[tree->nodeCount+7] = nextLeaf;
    tree->nodeCount += 8;
}


unsigned int assignParticleToTree(Particle *P, int ipart, Tree *tree, const double BOX[3]) {
    const unsigned int leaf = findLeafForPosition(P[ipart].Pos[0], P[ipart].Pos[1], P[ipart].Pos[2], tree, BOX);

    P[ipart].leafIndex = leaf;
    if (tree->particleCounts[leaf] == 0) {
        tree->firstParticle[leaf] = ipart;
    }
    ++ tree->particleCounts[leaf];

    return leaf;
}

unsigned int findLeafForPosition(const double x, const double y, const double z, const Tree *tree, const double BOX[3])
{
    for (unsigned int l = 0; l < tree->nodeCount; ++l) {
        if (nodeIsLeaf(tree, l) && coordInsideNode(x, y, z, tree->nodes[l], BOX)) {
            return l;
        }
    }
    fprintf(stderr, "No leaf found for position, aborting...");
    exit(1);
}

bool nodeIsLeaf(const Tree *tree, unsigned int leaf) {
    return tree->firstParticle[leaf] >= 0;
}

bool coordInsideNode(const double x, const double y, const double z, const Morton key, const double BOX[3]) {
    double sideLength[3];
    getNodeSize(sideLength, key, BOX);

    double kx, ky, kz;
    key2Coord(key, &kx, &ky, &kz, BOX);

    if (x < kx || x > kx + sideLength[0]) {
        return false;
    }
    if (y < ky || y > ky + sideLength[1]) {
        return false;
    }
    if (z < kz || z > kz + sideLength[2]) {
        return false;
    }

    return true;
}

void setParticleRangesInTree(Particle *particles, const int npart, Tree *tree) {
    unsigned int leafIndex = 0, oldLeafIndex = 0;

    for (int ipart = 0; ipart < npart; ++ipart) {
        leafIndex = particles[ipart].leafIndex;

        if (ipart == 0 || leafIndex != oldLeafIndex) {
            oldLeafIndex = leafIndex;

            tree->firstParticle[leafIndex] = ipart;
        }
    }
}

int findNGB(Particle *P, const int ipart, const double hsml, const Tree *tree, int *ngblist, const double BOX[3]) {
    const unsigned int topNodeIndex = getTopNodeIndex(P, ipart, hsml, tree, BOX);
    const int found = findNeighboursInNode(P, ipart, hsml, tree, ngblist, topNodeIndex);
    return found;
}

unsigned int getTopNodeIndex(const Particle *P, const int ipart, const double hsml, const Tree *tree,
                             const double *BOX) {
    unsigned int nodeIndex = P[ipart].leafIndex;
    Morton node = tree->nodes[nodeIndex];

    while (isNotRootNode(node) && !nodeSurroundsSphere(node, P[ipart].Pos, hsml, BOX)) {
        nodeIndex = getParentNode(nodeIndex, tree);
        node = tree->nodes[nodeIndex];
    }
    return nodeIndex;
}

bool isNotRootNode(Morton node) {
    return node.key != 0u;
}

bool nodeSurroundsSphere(Morton node, const double *center, const double radius, const double *BOX) {
    double nodeLowerCoords[3], nodeSideLength[3];
    nodeToBox(node, nodeLowerCoords, nodeSideLength, BOX);

    return sphereInsideBox(nodeLowerCoords, nodeSideLength, center, radius);
}

void nodeToBox(Morton node, double* lowerCoords, double* sideLength, const double BOX[3]) {
    key2Coord(node, &lowerCoords[0], &lowerCoords[1], &lowerCoords[2], BOX);
    for (int i = 0; i < 3; ++i) {
        sideLength[i] = BOX[i] / (1 << key2Depth(node));
    }
}

unsigned int getParentNode(unsigned int nodeIndex, const Tree* tree) {
    return tree->parentNodes[nodeIndex];
}

int findNeighboursInNode(Particle *P, const int ipart, const double hsml, const Tree *tree, int *ngblist,
                               unsigned int nodeIndex) {
    Morton node = tree->nodes[nodeIndex];
    unsigned int leafIndex = getFirstSubnodeInNode(nodeIndex, tree);
    while (! nodeIsLeaf(tree, leafIndex)) {
        leafIndex = getFirstSubnodeInNode(leafIndex, tree);
    }
    const unsigned int firstLeafIndex = leafIndex;
    Morton leaf;
    int found = 0;

    do {
        found += findNeighboursInLeaf(P, ipart, hsml, tree, ngblist, found, leafIndex);
        leafIndex = getNextLeaf(leafIndex, tree);

        // @todo really need this? Probably yes
        while (! nodeIsLeaf(tree, leafIndex)) {
            leafIndex = getFirstSubnodeInNode(leafIndex, tree);
        }

        leaf = tree->nodes[leafIndex];
    } while(nodeContainsLeaf(node, leaf) && leafIndex != firstLeafIndex);

    return found;
}

unsigned int getFirstSubnodeInNode(unsigned int nodeIndex, const Tree *tree) {
    return tree->nextNodes[nodeIndex];
}

int findNeighboursInLeaf(Particle *P, const int ipart, const double hsml, const Tree *tree, int *ngblist, int ingb,
                         unsigned int leafIndex) {
    assert(ingb < NGBMAX);

    int found = 0;
    //Assume there can be only 1 particle per leaf anyway
    const int possibleNeighbour = tree->firstParticle[leafIndex];

    if (possibleNeighbour == ipart) {
        // Do not count particle as its own neighbour
        return found;
    }

    double r2 = 0.0, d;
    for (int i = 0; i < 3; ++i) {
        d = P[ipart].Pos[i] - P[possibleNeighbour].Pos[i];
        r2 += d*d;
    }

    if (r2 <= hsml*hsml) {
        //Is a neighbour
        found = 1;
        ngblist[ingb] = possibleNeighbour;
    }

    return found;
}

unsigned int getNextLeaf(unsigned int leafIndex, const Tree *tree) {
    return tree->nextNodes[leafIndex];
}

// @todo Biggest bottleneck in ngb search (previously?)
bool nodeContainsLeaf(Morton node, Morton leaf) {
    Morton nextNode = translateToNextKey(node);

    // nextNode only makes sense if we do not span until the end of the box
    bool xRoot = isLastCoordInDimension(node.x, node.level);
    bool yRoot = isLastCoordInDimension(node.y, node.level);
    bool zRoot = isLastCoordInDimension(node.z, node.level);

    //if lower corner of leaf is smaller than higher corner of node bigger or equal than lower corner then it is inside
    if (!xRoot && (leaf.x < node.x || leaf.x >= nextNode.x)) {
        return false;
    }
    if (!yRoot && (leaf.y < node.y || leaf.y >= nextNode.y)) {
        return false;
    }
    if (!zRoot && (leaf.z < node.z || leaf.z >= nextNode.z)) {
        return false;
    }

    return true;
}

void freeTreeContents(Tree *tree) {
    free(tree->nextNodes);
    free(tree->parentNodes);
    free(tree->particleCounts);
    free(tree->firstParticle);
    free(tree->nodes);
}
