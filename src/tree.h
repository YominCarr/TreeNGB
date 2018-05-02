#ifndef TREENGB_TREE_H
#define TREENGB_TREE_H

#include <stdint.h>
#include "particle.h"

#define NGBMAX 1000
#define KEY uint64_t
#define COORD uint64_t
#define MAXDEPTH 10
#define MAXLEAVES (2 << (3 * MAXDEPTH))
#define MAXLEAFSIZE 5 //200^(1/3) rounded down

struct Tree
{
    KEY* leafs;
    int count;
};

struct Tree buildTree(struct Particle *P, const int npart, const double BOX[3]);

int findNGB(const int ipart, const double hsml, const struct Tree tree, int ngblist[NGBMAX]);


void buildTreeSerial(struct Particle *P, const int npart, struct Tree* tree, const double BOX[3]);


KEY coord2Key ( double x, double y, double z, const double BOX[3] );
void key2Coord ( KEY key, double *x, double *y, double *z, const double BOX[3] );

COORD translateCoordFromDouble(double c, double box);
double translateCoordToDouble(COORD c, double box);

void sortTree(struct Tree* tree);

#endif //TREENGB_TREE_H
