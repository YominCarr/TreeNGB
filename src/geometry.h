#ifndef TREENGB_GEOMETRY_H
#define TREENGB_GEOMETRY_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

bool sphereInsideBox(const double *boxLowerCoord, const double *sideLength,
                     const double *sphereCenter, const double radius);
bool coordInsideBox(const double *boxLowerCoord, const double *sideLength, const double *coord);
double getMaxRadiusForSphereInBox(const double *boxLowerCoord, const double *sideLength, const double *sphereCenter);

#ifdef __cplusplus
}
#endif

#endif //TREENGB_GEOMETRY_H
