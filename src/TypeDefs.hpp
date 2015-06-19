#ifndef TYPEDEFS_HDR
#define TYPEDEFS_HDR

#include <vector>

typedef unsigned int uint;
typedef const unsigned int cuint;
typedef const int cint;
typedef const float cfloat;
typedef const double cdouble;
typedef const bool cbool;

typedef std::vector<float> fvec;

enum ProjectMode {PM_PERSPECTIVE, PM_ORTHOGRAPHIC};
enum TraceMode {TM_SLAB, TM_GRID};
enum ScatterMode {SM_ISOTROPIC, SM_ANISOTROPIC};
enum IntegrateMode {IM_FORCED, IM_RANDOM};

#endif
