#ifndef TREENGB_PARTICLE_H
#define TREENGB_PARTICLE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include "morton.h"

typedef struct {
    double Pos[3];
    unsigned int leafIndex;
} Particle;

typedef struct Tree Tree;

Particle* createRandomParticles(const int N, const double BOX[3]);

void sortParticlesByKey(Particle* particles, const int N, Tree* tree);

bool compareParticles(Particle *particles, int i, int j, Tree* tree);
void swapParticles(Particle *particles, int i, int j);

#ifdef __cplusplus
}
#endif

#endif //TREENGB_PARTICLE_H
