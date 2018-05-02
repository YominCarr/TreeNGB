#include <malloc.h>
#include "particle.h"

struct Particle *createRandomParticles(int N) {
    struct Particle* P = malloc(N * sizeof(struct Particle));

    fprintf(stderr, "Implement 'createRandomParticles'\n");

    return P;
}
