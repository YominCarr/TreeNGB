#include <math.h>
#include "morton.h"

Morton coord2Key(const double x, const double y, const double z, const double BOX[3]) {
    Morton key;

    key.x = translateCoordFromDouble(x, BOX[0]);
    key.y = translateCoordFromDouble(y, BOX[1]);
    key.z = translateCoordFromDouble(z, BOX[2]);

    return key;
}

void key2Coord(const Morton key, double *x, double *y, double *z, const double BOX[3]) {
    *x = translateCoordToDouble(key.x, BOX[0]);
    *y = translateCoordToDouble(key.y, BOX[1]);
    *z = translateCoordToDouble(key.z, BOX[2]);
}

int key2Depth(const Morton key) {
    return key.level;
}

COORD translateCoordFromDouble(const double c, const double box) {
    return c / box * pow(2.0, MAXLEVEL);
}

double translateCoordToDouble(const COORD c, const double box) {
    return c * box / pow(2.0, MAXLEVEL);
}
