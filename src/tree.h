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
#define MAXDEPTH MAXLEVEL
#define MAXLEAVES (1 << 29) //More goes beyond uint in MAXLEAFS * sizeof(Morton))
#define MAXLEAFSIZE 1 //Needs to be one otherwise we can not find all previously added particles during node splitting

typedef struct {
    Morton *leaves;
    int *firstParticle; //Assume particles are sorted then leaf i contains particles
    int *particleCounts; //firstParticle[i] -> firstParticle[i] + particleCounts[i] - 1
    //@todo do I actually need to count if it's either 0 or 1?
    int leafCount;
} Tree;

Tree buildTree(Particle *P, const int npart, const double BOX[3]);
Tree initalizeTree();


void createRootNode(Tree *tree, const double *BOX);
void buildTreeSerial(Particle *P, const int npart, Tree *tree, const double BOX[3]);

bool treeHasSpaceForSplittingOnce(Tree* tree, int newDepth);

int findLeafForPosition(const double x, const double y, const double z, const Tree *tree, const double BOX[3]);
bool coordInsideNode(const double x, const double y, const double z, const Morton key, const double BOX[3]);
void getNodeSize(double* sideLength, const Morton key, const double BOX[3]);

void splitNode(Particle *P, const int l, Tree *tree, const double BOX[3]);
int assignParticleToTree(Particle *P, int ipart, Tree *tree, const double BOX[3]);

void setParticleRangesInTree(Particle *particles, const int npart, Tree *tree);
int getIndexOfLeaf(Morton leaf, Tree* tree);



int findNGB(Particle *P, const int ipart, const double hsml, const Tree *tree, int ngblist[NGBMAX], const double BOX[3]);

bool isNotRootNode(Morton node);
bool nodeBiggerThanSphere(Morton node, const double center[3], const double radius, const double BOX[3]);
void nodeToBox(Morton node, double* center, double* sideLength, const double BOX[3]);
Morton getParentNode(Morton node);

int findNeighboursInNode(Particle *P, const int ipart, const double hsml, const Tree *tree, int *ngblist,
                               Morton node);

Morton getFirstLeafInNode(Morton node);
Morton getNextLeaf(Morton leaf);
bool nodeContainsLeaf(Morton node, Morton leaf);

int findNeighboursInLeaf(Particle *P, const int ipart, const double hsml, const Tree *tree, int *ngblist, int ingb,
                         Morton leaf);


void freeTreeContents(Tree *tree);

#ifdef __cplusplus
}
#endif

#endif //TREENGB_TREE_H
