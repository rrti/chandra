#ifndef CONSTANTS_HDR
#define CONSTANTS_HDR

#include <cstdlib>
#include <cmath>

#include "TypeDefs.hpp"

// NOTE: we can't make epsilon much smaller, or
// numerical inaccuracy shows up in calculations
// like (50.0f - EPSILON) where the result then
// is indistinguishable from 50.0f (also depends
// on optimization flags)
static cfloat FRAND_MAX   = float(RAND_MAX);
static cfloat M_PI2       = M_PI * 2.0f;
static cfloat DEG2RAD     = M_PI / 180.0f;
static cfloat RAD2DEG     = 180.0f / M_PI;
static cfloat EPSILON     = 0.0001f;
static cuint  PXL_MAX_VAL = 255;

#endif
