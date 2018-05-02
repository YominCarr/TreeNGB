#include <malloc.h>
#include "particle.h"
#include "omp.h"
#include "proto.h"

struct Particle *createRandomParticles(const int N, const double BOX[3]) {
    struct Particle* P = malloc(N * sizeof(struct Particle));

#pragma omp parallel
    for (int i = 0; i < N; ++i) {
        P[i].Pos[0] = erand48 ( Omp.Seed ) * BOX[0];
        P[i].Pos[1] = erand48 ( Omp.Seed ) * BOX[1];
        P[i].Pos[2] = erand48 ( Omp.Seed ) * BOX[2];
    }

    return P;
}
