#include <stdio.h>
#include <float.h>
#include <math.h>
#include "geometry.h"

bool sphereInsideBox(const double *boxLowerCoord, const double *sideLength,
                     const double *sphereCenter, const double radius) {
    if (!coordInsideBox(boxLowerCoord, sideLength, sphereCenter)) {
        return false;
    }

    const double maxRadius = getMaxRadiusForSphereInBox(boxLowerCoord, sideLength, sphereCenter);
    return maxRadius >= radius;
}

bool coordInsideBox(const double *boxLowerCoord, const double *sideLength, const double *coord) {
    for (int i = 0; i < 3; ++i) {
        if (coord[i] < boxLowerCoord[i]) {
            return false;
        }
        if (coord[i] > boxLowerCoord[i] + sideLength[i]) {
            return false;
        }
    }
    return true;
}

double getMaxRadiusForSphereInBox(const double *boxLowerCoord, const double *sideLength, const double *sphereCenter) {
    double minDistanceToWall = DBL_MAX, left, right;
    for (int i = 0; i < 3; ++i) {
        left = fabs(sphereCenter[i] - boxLowerCoord[i]);
        right = fabs(sphereCenter[i] - (boxLowerCoord[i] + sideLength[i]));

        minDistanceToWall = fmin(minDistanceToWall, left);
        minDistanceToWall = fmin(minDistanceToWall, right);
    }

    return minDistanceToWall;
}
