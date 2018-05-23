#include <malloc.h>
#include "tree.h"
#include "particle.h"
#include "omp.h"
#include "proto.h"

struct OpenMP_infos Omp;

Particle *createRandomParticles(const int N, const double BOX[3]) {
    Particle* particles = malloc(N * sizeof(Particle));

#pragma omp parallel
    for (int i = 0; i < N; ++i) {
        particles[i].leafIndex = 0;

        particles[i].Pos[0] = erand48 ( Omp.Seed ) * BOX[0];
        particles[i].Pos[1] = erand48 ( Omp.Seed ) * BOX[1];
        particles[i].Pos[2] = erand48 ( Omp.Seed ) * BOX[2];
    }

    return particles;
}

// Bubblesort
void sortParticlesByKey(Particle *particles, const int N, Tree* tree) {
    int n = N;
    do{
        int newn = 1;
        for (int i=0; i<n-1; ++i){
            if (compareParticles(particles, i, i+1, tree)) {
                swapParticles(particles, i, i+1);
                newn = i+1;
            }
        }
        n = newn;
    } while (n > 1);
}

bool compareParticles(Particle *particles, int i, int j, Tree* tree) {
    return tree->nodes[particles[i].leafIndex].key > tree->nodes[particles[j].leafIndex].key;
}

void swapParticles(Particle *particles, int i, int j) {
    double p[3];
    unsigned int leafIndex;

    for (int k = 0; k < 3; ++k) {
        p[k] = particles[i].Pos[k];
    }
    leafIndex = particles[i].leafIndex;

    for (int k = 0; k < 3; ++k) {
        particles[i].Pos[k] = particles[j].Pos[k];
    }
    particles[i].leafIndex = particles[j].leafIndex;

    for (int k = 0; k < 3; ++k) {
        particles[j].Pos[k] = p[k];
    }
    particles[j].leafIndex = leafIndex;
}
