#ifndef GRID_HDR
#define GRID_HDR

#include <vector>

#include "TypeDefs.hpp"
#include "vec3.hpp"

struct RNG;
struct Photon;

struct Grid {
	Grid(cuint _xres, cuint _yres, cuint _zres, cfloat _xsize, cfloat _ysize, cfloat _zsize, cfloat _den, cfloat _opa);
	void Init();

	float GetMaxOpticalDepth(const Photon& p) const;
	void SetForcedFirstScatterPos(Photon& p, const RNG& rng) const;
	bool Integrate(Photon& p, cfloat tau, IntegrateMode im) const;

	// these calculate maximum distance photon
	// can travel in each "cardinal" direction
	float GetMaxXDistance(const Photon& pho, const vec3& pcur) const;
	float GetMaxYDistance(const Photon& pho, const vec3& pcur) const;
	float GetMaxZDistance(const Photon& pho, const vec3& pcur) const;

	// these map 3D POSITION <0, 0, 0> to INDICES
	// <xres/ 2 - 1, yres / 2 - 1, zres / 2 - 1>
	// (as if the position <x = 0, y = 0, z = 0>
	// was <x + hsize.x, y + hsize.y, z + hsize.z>)
	//
	// cell boundaries lie in [0.0, size], but they represent
	// 3D coordinates in [-hsize, +hsize] (ie. "world-space")
	int GetXCellIdx(cfloat x) const { return int(xres * (x + hsize.x) / xsize); }
	int GetYCellIdx(cfloat y) const { return int(yres * (y + hsize.y) / ysize); }
	int GetZCellIdx(cfloat z) const { return int(zres * (z + hsize.z) / zsize); }

	bool IndicesInBounds(cint i, cint j, cint k) const;
	bool RPosInBounds(const vec3& p) const;
	bool GPosInBounds(const vec3& p) const;
	void ClampRPos(vec3& p) const;
	void ClampGPos(vec3& p) const;

	uint xres, yres, zres;						// number of cells in each dimension (resolution)
	float hxsize, hysize, hzsize;				// grid half-size (coor. span) in each dimension
	float  xsize,  ysize,  zsize;				// grid full-size (coor. span) in each dimension
	fvec xfaces, yfaces, zfaces;				// cell faces (boundaries) along each axis
	std::vector<std::vector<fvec > > rhokappa;	// 3D cell density array

	vec3 hsize;									// (hxsize, hysize, hzsize)
	vec3 fsize;									// ( xsize,  ysize,  zsize)

	float density;
	float opacity;
};

#endif
