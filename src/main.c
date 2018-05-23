#include <stdio.h>
#include <malloc.h>
#include <omp.h>
#include <stdlib.h>
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

    const int NPART = 1000;
    double BOX[3] = {1.0, 1.0, 1.0};
    printf("Creating %d random particles inside box of %g x %g x %g\n\n", NPART, BOX[0], BOX[1], BOX[2]);

    Particle *P = createRandomParticles(NPART, BOX);

    printf("Build tree\n\n");
    Tree tree = buildTree(P, NPART, BOX);

    printf("Search for neighbours\n\n");
    int* ngblist = calloc(NGBMAX, sizeof(int));
    int found = findNGB(P, 0, 0.1, &tree, ngblist, BOX);
    printf("Found %d neighbours of particle 0 in distance 0.1\n\n", found);

    freeTreeContents(&tree);
    free(P);

    return 0;
}

