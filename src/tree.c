#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tree.h"
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

    tree.leafs = calloc(MAXLEAFS, sizeof(Morton));
    tree.firstParticle = calloc(MAXLEAFS, sizeof(int));
    tree.particleCounts = calloc(MAXLEAFS, sizeof(int));

    tree.leafCount = 0;

    return tree;
}

int findNGB(const int ipart, const double hsml, const Tree tree, int *ngblist) {
    fprintf(stderr, "Implement 'findNGB'\n");
    return 0;
}

void createRootNode(Tree *tree, const double *BOX) {
    tree->leafCount = 1;

    const Morton node = coord2Key(0.0, 0.0, 0.0, BOX);
    tree->leafs[0] = node;

}

void buildTreeSerial(Particle *P, const int npart, Tree *tree, const double *BOX) {
    for (int ipart = 0; ipart < npart; ++ipart) {
        int leaf = assignParticleToTree(P, ipart, tree, BOX);

        while (tree->particleCounts[leaf] > MAXLEAFSIZE) {
            if (key2Depth(tree->leafs[leaf]) == MAXDEPTH) {
                fprintf(stderr, "Trying to refine tree beyond capacity, aborting...");
                exit(1);
            }

            splitNode(P, leaf, tree, BOX);
            leaf = assignParticleToTree(P, ipart, tree, BOX);
        }
    }

    sortParticlesByKey(P, npart);
    //After resorting particles we need to readjust the tree
    // @todo maybe absorb this directly into the resorting?
    setParticleRangesInTree(P, npart, tree);
}

void splitNode(Particle *P, const int l, Tree *tree, const double BOX[3]) {
    const Morton parent = tree->leafs[l];
    double pX, pY, pZ; //Assume these are the corner with the smallest coord
    key2Coord(parent, &pX, &pY, &pZ, BOX);
    const int parentLevel = key2Depth(parent);

    //Particle to redistribute
    const int epart = tree->firstParticle[l];

    //Clean node spot which will be reused
    tree->particleCounts[l] = 0;

    double newSize[3];
    for (int i = 0; i < 3; ++i) {
        newSize[i] = BOX[i] / (2 << parentLevel);
    }

    //Create 8 new nodes, save one in l and the others in Tree.leafCount (+0,1,2,...)
    int save[8] = {l, tree->leafCount, tree->leafCount+1, tree->leafCount+2, tree->leafCount+3,
                   tree->leafCount+4, tree->leafCount+5, tree->leafCount+6};
    for (int i = 0, c = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k, ++c) {
                const double x = pX + i * newSize[0];
                const double y = pY + j * newSize[1];
                const double z = pZ + k * newSize[2];
                Morton node = coord2Key(x, y, z, BOX);
                node.level = parentLevel + 1;

                const int s = save[c];
                tree->leafs[s] = node;

                if (coordInsideNode(P[epart].Pos[0], P[epart].Pos[1], P[epart].Pos[2], node, BOX)) {
                    P[epart].leaf = node;
                    tree->particleCounts[s] = 1;
                }
            }
        }
    }
    tree->leafCount += 7;
}


int assignParticleToTree(Particle *P, int ipart, Tree *tree, const double BOX[3]) {
    const int leaf = findLeafForPosition(P[ipart].Pos[0], P[ipart].Pos[1], P[ipart].Pos[2], tree, BOX);

    P[ipart].leaf = tree->leafs[leaf];
    ++ tree->particleCounts[leaf];

    return leaf;
}

int findLeafForPosition(const double x, const double y, const double z, const Tree *tree, const double BOX[3])
{
    for (int l = 0; l < tree->leafCount; ++l) {
        if (coordInsideNode(x, y, z, tree->leafs[l], BOX)) {
            return l;
        }
    }
    fprintf(stderr, "No leaf found for position, aborting...");
    exit(1);
}

bool coordInsideNode(const double x, const double y, const double z, const Morton key, const double *BOX) {
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

void getNodeSize(double* sideLength, const Morton key, const double BOX[3]) {
    for (int i = 0; i < 3; ++i) {
        sideLength[i] = BOX[i] / pow(2.0, key2Depth(key));
    }
}

void setParticleRangesInTree(Particle *particles, const int npart, Tree *tree) {
    Morton leaf, oldLeaf;
    oldLeaf.key = 0;
    int leafIndex = 0;

    for (int ipart = 0; ipart < npart; ++ipart) {
        leaf = particles[ipart].leaf;

        if (ipart == 0 || leaf.key != oldLeaf.key) {
            leafIndex = getIndexOfLeaf(leaf, tree);
            oldLeaf = leaf;

            tree->firstParticle[leafIndex] = ipart;
        }
    }
}

int getIndexOfLeaf(Morton leaf, Tree* tree) {
    for (int i = 0; i < tree->leafCount; ++i) {
        if (tree->leafs[i].key == leaf.key) {
            return i;
        }
    }
    fprintf(stderr, "Leaf %lu unknown, aborting...", leaf.key);
    exit(1);
}

void freeTreeContents(Tree *tree) {
    free(tree->particleCounts);
    free(tree->firstParticle);
    free(tree->leafs);
}
