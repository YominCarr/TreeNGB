#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include "tree.h"
#include "external/morton_utils.h"

struct Tree buildTree(struct Particle *P, const int npart, const double BOX[3]) {
    fprintf(stderr, "Implement a parallel tree build\n");

    struct Tree tree;
    tree.leafs = malloc(MAXLEAVES * sizeof(KEY));
    tree.particleCounts = malloc(MAXLEAVES * sizeof(int));
    tree.leafCount = 0;

    buildTreeSerial(P, npart, &tree, BOX);
    sortTree(&tree);

    return tree;
}

int findNGB(const int ipart, const double hsml, const struct Tree tree, int *ngblist) {
    fprintf(stderr, "Implement 'findNGB'\n");
    return 0;
}

void buildTreeSerial(struct Particle *P, const int npart, struct Tree *tree, const double *BOX) {
    for (int ipart = 0; ipart < npart; ++ipart) {
        const int l = findLeafForPosition(P[ipart].Pos[0], P[ipart].Pos[1], P[ipart].Pos[2], tree, BOX);

        P[ipart].treeIndex = l;
        ++ tree->particleCounts[l];
        if (tree->particleCounts[l] > MAXLEAFSIZE && key2Depth(tree->leafs[l] < MAXDEPTH)) {
            splitNode(l, tree, BOX);
        }
    }
}

void splitNode(const int l, struct Tree *tree, const double BOX[3]) {
    fprintf(stderr, "Implement 'splitNode'\n");
    const KEY parent = tree->leafs[l];
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
                const KEY node = coord2Key(x, y, z, BOX);
                const int s = save[c];
                tree->leafs[s] = node;
                //@todo assign particles to 8 new nodes and set tree->particleCounts[s]
                //@todo how do I get all particles inside node in a quick way? Need that for ngb search also
            }
        }
    }
    tree->leafs += 7;
}

int findLeafForPosition(const double x, const double y, const double z, const struct Tree *tree, const double BOX[3])
{
    for (int l = 0; l < tree->leafCount; ++l) {
        if (coordInsideNode(x, y, z, tree->leafs[l], BOX)) {
            return l;
        }
    }
    fprintf(stderr, "No leaf found for position, aborting...");
    exit(1);
}

bool coordInsideNode(const double x, const double y, const double z, const KEY key, const double *BOX) {
    fprintf(stderr, "Implement 'coordInsideNode'\n");
    return 0;
}

KEY coord2Key(const double x, const double y, const double z, const double BOX[3]) {
    const COORD xT = translateCoordFromDouble(x, BOX[0]);
    const COORD yT = translateCoordFromDouble(y, BOX[1]);
    const COORD zT = translateCoordFromDouble(z, BOX[2]);

    return coord2morton3D_64(xT, yT, zT);
}

void key2Coord(const KEY key, double *x, double *y, double *z, const double BOX[3]) {
    COORD xT, yT, zT;
    morton2coord3D_64(key, &xT, &yT, &zT);

    *x = translateCoordToDouble(xT, BOX[0]);
    *y = translateCoordToDouble(yT, BOX[1]);
    *z = translateCoordToDouble(zT, BOX[2]);
}

int key2Depth(const KEY key) {
    fprintf(stderr, "Implement 'key2Depth'\n");
    return -1;
}

COORD translateCoordFromDouble(const double c, const double box) {
    fprintf(stderr, "Implement 'translateCoordFromDouble'\n");
    return 0;
}

double translateCoordToDouble(const COORD c, const double box) {
    fprintf(stderr, "Implement 'translateCoordToDouble'\n");
    return 0;
}

void sortTree(struct Tree *tree) {
    fprintf(stderr, "Implement 'sortTree'\n");
}
