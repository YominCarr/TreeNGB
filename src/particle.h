#ifndef TREENGB_PARTICLE_H
#define TREENGB_PARTICLE_H

#include <stdbool.h>
#include "morton.h"

typedef struct {
    double Pos[3];
    Morton leaf;
} Particle;

Particle* createRandomParticles(const int N, const double BOX[3]);

void sortParticlesByKey(Particle* particles, const int N);

bool compareParticles(Particle *particles, int i, int j);
void swapParticles(Particle *particles, int i, int j);

#endif //TREENGB_PARTICLE_H
