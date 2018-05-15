#include <malloc.h>
#include "particle.h"
#include "omp.h"
#include "proto.h"

struct OpenMP_infos Omp;

Particle *createRandomParticles(const int N, const double BOX[3]) {
    Particle* particles = malloc(N * sizeof(Particle));

#pragma omp parallel
    for (int i = 0; i < N; ++i) {
        particles[i].leaf.key = 0;

        particles[i].Pos[0] = erand48 ( Omp.Seed ) * BOX[0];
        particles[i].Pos[1] = erand48 ( Omp.Seed ) * BOX[1];
        particles[i].Pos[2] = erand48 ( Omp.Seed ) * BOX[2];
    }

    return particles;
}

// Bubblesort
void sortParticlesByKey(Particle *particles, const int N) {
    int n = N;
    do{
        int newn = 1;
        for (int i=0; i<n-1; ++i){
            if (compareParticles(particles, i, i+1)) {
                swapParticles(particles, i, i+1);
                newn = i+1;
            }
        }
        n = newn;
    } while (n > 1);
}

bool compareParticles(Particle *particles, int i, int j) {
    return particles[i].leaf.key > particles[j].leaf.key;
}

void swapParticles(Particle *particles, int i, int j) {
    fprintf(stderr, "Implement 'swapParticles'\n");
}
