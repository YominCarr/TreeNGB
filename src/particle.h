#ifndef TREENGB_PARTICLE_H
#define TREENGB_PARTICLE_H

struct Particle {
    double Pos[3];
};

struct Particle* createRandomParticles(const int N, const double BOX[3]);

#endif //TREENGB_PARTICLE_H
