#ifndef TREENGB_MORTON_H
#define TREENGB_MORTON_H

#include <stdint.h>

typedef uint32_t COORD;

typedef union {
    uint64_t key;
    struct { //Partitioning to allow the deepest level: max 19
        int x:19, y:19, z:19;
        int level:6;
        int assigned:1;
    };
} Morton;

#define MAXLEVEL 19

Morton coord2Key(const double x, const double y, const double z, const double BOX[3]);
void key2Coord(const Morton key, double *x, double *y, double *z, const double BOX[3]);
int key2Depth(const Morton key);

COORD translateCoordFromDouble(const double c, const double box);
double translateCoordToDouble(const COORD c, const double box);

#endif //TREENGB_MORTON_H
