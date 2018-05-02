#include <stdio.h>
#include <malloc.h>
#include "tree.h"
#include "external/morton_utils.h"

struct Tree buildTree(struct Particle *P, const int npart, const double BOX[3]) {
    fprintf(stderr, "Implement a parallel tree build\n");

    struct Tree tree;
    tree.leafs = malloc(MAXLEAVES * sizeof(KEY));
    tree.count = 0;

    buildTreeSerial(P, npart, &tree, BOX);
    sortTree(&tree);

    return tree;
}

int findNGB(const int ipart, const double hsml, const struct Tree tree, int *ngblist) {
    fprintf(stderr, "Implement 'findNGB'\n");
    return 0;
}

void buildTreeSerial(struct Particle *P, const int npart, struct Tree *tree, const double *BOX) {
    fprintf(stderr, "Implement 'buildTreeSerial'\n");
    for (int ipart = 0; ipart < npart; ++ipart) {

    }
}

KEY coord2Key(double x, double y, double z, const double BOX[3]) {
    const COORD xT = translateCoordFromDouble(x, BOX[0]);
    const COORD yT = translateCoordFromDouble(y, BOX[1]);
    const COORD zT = translateCoordFromDouble(z, BOX[2]);

    return coord2morton3D_64(xT, yT, zT);
}

void key2Coord(KEY key, double *x, double *y, double *z, const double BOX[3]) {
    COORD xT, yT, zT;
    morton2coord3D_64(key, &xT, &yT, &zT);

    *x = translateCoordToDouble(xT, BOX[0]);
    *y = translateCoordToDouble(yT, BOX[1]);
    *z = translateCoordToDouble(zT, BOX[2]);
}

COORD translateCoordFromDouble(double c, double box) {
    fprintf(stderr, "Implement 'translateCoordFromDouble'\n");
    return 0;
}

double translateCoordToDouble(COORD c, double box) {
    fprintf(stderr, "Implement 'translateCoordToDouble'\n");
    return 0;
}

void sortTree(struct Tree *tree) {
    fprintf(stderr, "Implement 'sortTree'\n");
}
