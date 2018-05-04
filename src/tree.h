#ifndef TREENGB_TREE_H
#define TREENGB_TREE_H

#include <stdint.h>
#include <stdbool.h>
#include "morton.h"
#include "particle.h"

#define NGBMAX 1000
#define MAXDEPTH 10 //More goes beyond integer in MAXLEAVES
#define MAXLEAVES (2 << (3 * MAXDEPTH))
#define MAXLEAFSIZE 5 //200^(1/3) rounded down

typedef struct {
    Morton *leafs;
    int *firstParticle; //Assume particles are sorted then leaf i contains particles
    int *particleCounts; //firstParticle[i] -> firstParticle[i] + particleCounts[i] - 1
    int leafCount;
} Tree;

Tree buildTree(Particle *P, const int npart, const double BOX[3]);

int findNGB(const int ipart, const double hsml, const Tree tree, int ngblist[NGBMAX]);


void buildTreeSerial(Particle *P, const int npart, Tree *tree, const double BOX[3]);

int findLeafForPosition(const double x, const double y, const double z, const Tree *tree, const double BOX[3]);
bool coordInsideNode(const double x, const double y, const double z, const Morton key, const double BOX[3]);

void splitNode(const int l, Tree *tree, const double BOX[3]);

void setParticleRangesInTree(Particle *particles, const int npart, Tree *tree);
int getIndexOfLeaf(Morton leaf, Tree* tree);

void freeTreeContents(Tree *tree);

#endif //TREENGB_TREE_H
