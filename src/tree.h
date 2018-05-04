#ifndef TREENGB_TREE_H
#define TREENGB_TREE_H

#include <stdint.h>
#include <stdbool.h>
#include "particle.h"
#include "morton.h"

#define NGBMAX 1000
#define MAXDEPTH 10 //More goes beyond integer in MAXLEAVES
#define MAXLEAVES (2 << (3 * MAXDEPTH))
#define MAXLEAFSIZE 5 //200^(1/3) rounded down

typedef struct {
    Morton *leafs;
    int *particleCounts;
    int leafCount;
} Tree;

Tree buildTree(struct Particle *P, const int npart, const double BOX[3]);

int findNGB(const int ipart, const double hsml, const Tree tree, int ngblist[NGBMAX]);


void buildTreeSerial(struct Particle *P, const int npart, Tree *tree, const double BOX[3]);

int findLeafForPosition(const double x, const double y, const double z, const Tree *tree, const double BOX[3]);
bool coordInsideNode(const double x, const double y, const double z, const Morton key, const double BOX[3]);

void splitNode(const int l, Tree *tree, const double BOX[3]);

void sortTree(Tree *tree);

#endif //TREENGB_TREE_H
