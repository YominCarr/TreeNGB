#ifndef TREENGB_TREE_H
#define TREENGB_TREE_H

#include <stdint.h>
#include <stdbool.h>
#include "particle.h"

#define NGBMAX 1000
#define KEY uint64_t
#define COORD uint64_t
#define MAXDEPTH 10
#define MAXLEAVES (2 << (3 * MAXDEPTH))
#define MAXLEAFSIZE 5 //200^(1/3) rounded down

typedef struct {
    KEY *leafs;
    int *particleCounts;
    int leafCount;
} Tree;

Tree buildTree(struct Particle *P, const int npart, const double BOX[3]);

int findNGB(const int ipart, const double hsml, const Tree tree, int ngblist[NGBMAX]);


void buildTreeSerial(struct Particle *P, const int npart, Tree *tree, const double BOX[3]);

int findLeafForPosition(const double x, const double y, const double z, const Tree *tree, const double BOX[3]);
bool coordInsideNode(const double x, const double y, const double z, const KEY key, const double BOX[3]);

void splitNode(const int l, Tree *tree, const double BOX[3]);

KEY coord2Key(const double x, const double y, const double z, const double BOX[3]);
void key2Coord(const KEY key, double *x, double *y, double *z, const double BOX[3]);
int key2Depth(const KEY key);

COORD translateCoordFromDouble(const double c, const double box);
double translateCoordToDouble(const COORD c, const double box);

void sortTree(Tree *tree);

#endif //TREENGB_TREE_H
