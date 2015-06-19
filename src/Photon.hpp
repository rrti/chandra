#ifndef PHOTON_HDR
#define PHOTON_HDR

#include "TypeDefs.hpp"
#include "vec3.hpp"

struct Dust;
struct Slab;
struct Grid;
struct Image;
struct Observer;
struct RNG;

struct Photon {
	Photon(const RNG& rng, const vec3& p, const Grid& g);

	void Scatter(const RNG& rng, const Dust& dust, ScatterMode sm);
	void IsotropicScatter(const RNG& rng, const Dust&);
	void AnisotropicScatter(const RNG& rng, const Dust& dust);
	void AnisotropicScatterBW(const RNG& rng, const Dust& dust);

	void UpdateAngles(cfloat t, cfloat p);
	void UpdateDir();
	void UpdatePos(cfloat d);
	void UpdateCellIndices(const Grid& g);
	std::string str();

	bool PerspProject(const Observer& o, Image* img) const;
	bool OrthoProject(const Observer& o, Image* img) const;

	vec3 pos;				// position in 3D Cartesian coors
	vec3 dir;				// direction in 3D Cartesian coors
	vec4 pol;				// Stokes polarization state-vector (x=I, y=Q, z=U, w=V)
	float the, phi;			// zenith, azimuth angle in radians
	float cthe, sthe;		// {cos, s}ine of zenith angle in spherical coors
	float cphi, sphi;		// {cos, s}ine of azimuth angle in spherical coors

	uint numScatEvents;
	bool absorbed;

	int xcell;
	int ycell;
	int zcell;
};

#endif
