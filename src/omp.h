#ifndef TREENGB_OMP_H
#define TREENGB_OMP_H

extern struct OpenMP_infos {
    int NThreads;          // Number of openMP threads
    int ThreadID;          // Thread ID of this thread
    unsigned short Seed[3]; // random number seed: erand48(Omp.Seed)
} Omp;
#pragma omp threadprivate(Omp)

#endif //TREENGB_OMP_H
