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

struct Tree {
    Morton *nodes;
    int *firstParticle;
    int *particleCounts; // if negative the node is not a leaf
    unsigned int *parentNodes;
    //@todo do I actually need to count if it's either 0 or 1?
    int nodeCount;
};

typedef struct Tree Tree;

Tree buildTree(Particle *P, const int npart, const double BOX[3]);
Tree initalizeTree();


void createRootNode(Tree *tree, const double *BOX);
void buildTreeSerial(Particle *P, const int npart, Tree *tree, const double BOX[3]);

bool treeHasSpaceForSplittingOnce(Tree* tree, int newDepth);

unsigned int findLeafForPosition(const double x, const double y, const double z, const Tree *tree, const double BOX[3]);
bool nodeIsLeaf(const Tree *tree, unsigned int leaf);
bool coordInsideNode(const double x, const double y, const double z, const Morton key, const double BOX[3]);
void getNodeSize(double* sideLength, const Morton key, const double BOX[3]);

void splitNode(Particle *P, const unsigned int l, Tree *tree, const double BOX[3]);
unsigned int assignParticleToTree(Particle *P, int ipart, Tree *tree, const double BOX[3]);

void setParticleRangesInTree(Particle *particles, const int npart, Tree *tree);



int findNGB(Particle *P, const int ipart, const double hsml, const Tree *tree, int ngblist[NGBMAX], const double BOX[3]);

bool isNotRootNode(Morton node);
bool nodeSurroundsSphere(Morton node, const double *center, const double radius, const double *BOX);
void nodeToBox(Morton node, double* lowerCoords, double* sideLength, const double BOX[3]);
unsigned int getParentNode(unsigned int nodeIndex, const Tree *tree);

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
