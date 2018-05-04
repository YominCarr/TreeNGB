#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include "tree.h"
#include "morton.h"

Tree buildTree(Particle *P, const int npart, const double BOX[3]) {
    fprintf(stderr, "Implement a parallel tree build\n");

    Tree tree;
    tree.leafs = malloc(MAXLEAVES * sizeof(Morton));
    tree.firstParticle = malloc(MAXLEAVES * sizeof(int));
    tree.particleCounts = malloc(MAXLEAVES * sizeof(int));
    tree.leafCount = 0;

    buildTreeSerial(P, npart, &tree, BOX);

    return tree;
}

int findNGB(const int ipart, const double hsml, const Tree tree, int *ngblist) {
    fprintf(stderr, "Implement 'findNGB'\n");
    return 0;
}

void buildTreeSerial(Particle *P, const int npart, Tree *tree, const double *BOX) {
    for (int ipart = 0; ipart < npart; ++ipart) {
        const int l = findLeafForPosition(P[ipart].Pos[0], P[ipart].Pos[1], P[ipart].Pos[2], tree, BOX);

        P[ipart].leaf = tree->leafs[l];
        ++ tree->particleCounts[l];

        if (tree->particleCounts[l] > MAXLEAFSIZE && key2Depth(tree->leafs[l]) < MAXDEPTH) {
            splitNode(l, tree, BOX);
        }
    }

    sortParticlesByKey(P, npart);
    setParticleRangesInTree(P, npart, tree);
}

void splitNode(const int l, Tree *tree, const double BOX[3]) {
    fprintf(stderr, "Implement 'splitNode'\n");
    const Morton parent = tree->leafs[l];
    double pX, pY, pZ; //Assume these are the corner with the smallest coord
    key2Coord(parent, &pX, &pY, &pZ, BOX);
    const int parentLevel = key2Depth(parent);

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
                const Morton node = coord2Key(x, y, z, BOX);
                const int s = save[c];
                tree->leafs[s] = node;
                //@todo assign particles to 8 new nodes and set tree->particleCounts[s]
                //@todo use firstParticle field assuming particles are sorted (actually this is done later)
                //@todo how do I know which particles have to be put in which subnode?
            }
        }
    }
    tree->leafs += 7;
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
    fprintf(stderr, "Implement 'coordInsideNode'\n");
    return 0;
}

void setParticleRangesInTree(Particle *particles, const int npart, Tree *tree) {
    //@todo also set particle count here? maybe already done
    Morton leaf, oldLeaf;
    oldLeaf.key = 0;
    int leafIndex;

    for (int ipart = 0; ipart < npart; ++ipart) {
        leaf = particles[ipart].leaf;

        if (ipart == 0 || leaf.key != oldLeaf.key) {
            leafIndex = getIndexOfLeaf(leaf, tree);
            oldLeaf = leaf;

            tree->firstParticle[leafIndex] = ipart;
            tree->particleCounts[leafIndex] = 1;
        } else {
            ++ tree->particleCounts[leafIndex];
        }
    }
}

int getIndexOfLeaf(Morton leaf, Tree* tree) {
    fprintf(stderr, "Implement 'getIndexOfLeaf'\n");
    return 0;
}

void freeTreeContents(Tree *tree) {
    free(tree->particleCounts);
    free(tree->firstParticle);
    free(tree->leafs);
}
