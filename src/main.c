#include <stdio.h>
#include <malloc.h>
#include <omp.h>
#include "particle.h"
#include "tree.h"
#include "omp.h"
#include "proto.h"

#pragma omp threadprivate(Omp)
struct OpenMP_infos Omp = {0};

int main(int argc, char *argv[]) {
    printf("Welcome to TreeNGB...\n\n");

#pragma omp parallel
    {

        Omp.ThreadID = omp_get_thread_num();
        Omp.NThreads = omp_get_num_threads();
        Omp.Seed[2] = 14041981 * (Omp.ThreadID + 1);
        erand48(Omp.Seed); // remove leading 0
    }

    const int N = 1000;
    double BOX[3] = {1.0, 1.0, 1.0};
    printf("Creating %d random particles inside box of %g x %g x %g\n\n", N, BOX[0], BOX[1], BOX[2]);

    struct Particle *P = createRandomParticles(N, BOX);

    printf("Build tree\n\n");
    buildTree();

    free(P);

    return 0;
}

