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
    return c / box * (1 << MAXLEVEL);
}

double translateCoordToDouble(const COORD c, const double box) {
    return c * box / (1 << MAXLEVEL);
}

// @todo currently the easiest spot to save time in tree build and ngb search if we save the sizes instead of on demand
void getNodeSize(double* sideLength, const Morton key, const double BOX[3]) {
    const int depth = key2Depth(key);
    const double fac = 1.0 / (1 << depth);

    for (int i = 0; i < 3; ++i) {
        sideLength[i] = BOX[i] * fac;
    }
}

Morton translateToNextKey(const Morton key) {
    Morton newKey;

    newKey.x = translateToNextCoord(key.x, key.level);
    newKey.y = translateToNextCoord(key.y, key.level);
    newKey.z = translateToNextCoord(key.z, key.level);
    newKey.level = key.level;

    return newKey;
}

COORD translateToNextCoord(const COORD c, const unsigned int level) {
    return c + (1 << (MAXLEVEL - level));
}

bool isLastCoordInDimension(const COORD c, const unsigned int level) {
    const unsigned int max = (1 << MAXLEVEL);
    const unsigned int difference = (1 << (MAXLEVEL - level));
    return c == max - difference;
}
