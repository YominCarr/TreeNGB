#include <stdio.h>
#include "geometry.h"

bool sphereInsideBox(const double *boxCenter, const double *sideLength,
                     const double *sphereCenter, const double radius) {
    if (! sphereCenterInsideBox(boxCenter, sideLength, sphereCenter)) {
        return false;
    }

    const double maxRadius = getMaxRadiusForSphereInBox(boxCenter, sideLength, sphereCenter);
    return maxRadius >= radius;
}

bool sphereCenterInsideBox(const double *boxCenter, const double *sideLength, const double *sphereCenter) {
    fprintf(stderr, "Implement sphereCenterInsideBox!\n");
    return 0;
}

double getMaxRadiusForSphereInBox(const double *boxCenter, const double *sideLength, const double *sphereCenter) {
    fprintf(stderr, "Implement getMaxRadiusForSphereInBox!\n");
    return 0;
}
