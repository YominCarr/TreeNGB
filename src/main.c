#include <stdio.h>
#include <malloc.h>
#include "particle.h"
#include "tree.h"

int main ( int argc, char *argv[] )
{
    printf("Welcome to TreeNGB...\n\n");

    const int N = 1000;
    double BOX[3] = {1.0, 1.0, 1.0};
    printf("Creating %d random particles inside box of %g x %g x %g\n\n", N, BOX[0], BOX[1], BOX[2]);

    struct Particle* P = createRandomParticles(N);

    printf("Build tree\n\n");
    buildTree();

    free(P);

    return 0;
}

