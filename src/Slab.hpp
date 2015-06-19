#ifndef SLAB_HDR
#define SLAB_HDR

#include "TypeDefs.hpp"
#include "vec3.hpp"

class CTDFStackParser;
struct RNG;
struct ThreadData;
struct Photon;

struct Slab {
	Slab(cfloat _zmin, cfloat _zmax, cfloat _den, cfloat _opa, cuint _numLvls);

	void UpdateIntensityMoments(ThreadData* td, const vec3& opos, const Photon& p);

	inline bool PosInBounds(const vec3& p) const {
		return (p.z >= zmin && p.z <= zmax);
	}

	float GetRandomL(const RNG& rng) const;

	// depth is (zmax - zmin)
	float zmin;
	float zmax;
	float zRange;
	float density;
	float opacity;
	float tauMax;

	uint numLevels;
};

#endif
