#ifndef TREENGB_PARTICLE_H
#define TREENGB_PARTICLE_H

struct Particle {
    double Pos[3];
};

struct Particle* createRandomParticles(int N);

#endif //TREENGB_PARTICLE_H
