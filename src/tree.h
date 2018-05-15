#ifndef TREENGB_TREE_H
#define TREENGB_TREE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stdbool.h>
#include "morton.h"
#include "particle.h"

#define NGBMAX 1000
#define MAXDEPTH 9 //More goes beyond uint in MAXLEAVES * sizeof(Morton))
#define MAXLEAFS (1 << (3 * MAXDEPTH))
// @todo actually the tree may be deeper at some points, since we are filling up from the left
#define MAXLEAFSIZE 1 //Needs to be one otherwise we can not find all previously added particles during node splitting
//@todo looks like maybe morton keys can be even only 32 instead 64 of size

typedef struct {
    Morton *leafs;
    int *firstParticle; //Assume particles are sorted then leaf i contains particles
    int *particleCounts; //firstParticle[i] -> firstParticle[i] + particleCounts[i] - 1
    //@todo do I actually need to count if it's either 0 or 1?
    int leafCount;
} Tree;

Tree buildTree(Particle *P, const int npart, const double BOX[3]);
Tree initalizeTree();

int findNGB(const int ipart, const double hsml, const Tree tree, int ngblist[NGBMAX]);


void createRootNode(Tree *tree, const double *BOX);
void buildTreeSerial(Particle *P, const int npart, Tree *tree, const double BOX[3]);

int findLeafForPosition(const double x, const double y, const double z, const Tree *tree, const double BOX[3]);
bool coordInsideNode(const double x, const double y, const double z, const Morton key, const double BOX[3]);

void splitNode(Particle *P, const int l, Tree *tree, const double BOX[3]);
int assignParticleToTree(Particle *P, int ipart, Tree *tree, const double BOX[3]);

void setParticleRangesInTree(Particle *particles, const int npart, Tree *tree);
int getIndexOfLeaf(Morton leaf, Tree* tree);

void freeTreeContents(Tree *tree);

#ifdef __cplusplus
}
#endif

#endif //TREENGB_TREE_H
