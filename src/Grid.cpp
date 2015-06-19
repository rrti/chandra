#include <fstream>
#include <vector>
#include <cassert>

#include "TypeDefs.hpp"
#include "RNG.hpp"
#include "Grid.hpp"
#include "Photon.hpp"

// NOTE: if the number of cells is
// even in a given dimension, then
// no clear "middle" cell defined
//
// NOTE: grid origin (corner) is at (0, 0, 0), but
// every photon starting at (0, 0, 0) has indices
// (xres / 2 - 1, yres / 2 - 1, zres / 2 - 1) and
// is traced as if the middle cell is its starting
// point
//
// grid itself lives in the positive (+x, +y, +z)
// quadrant from (0, 0, 0) to (xsize, ysize, zsize)
// (gnuplot the volume .dat file)
//
// x ranges from 0 to (xres * xsize) / xres: (0, xsize)
// y ranges from 0 to (yres * ysize) / yres: (0, ysize)
// z ranges from 0 to (zres * zsize) / zres: (0, zsize)
//
// WARNING: when grid is extremely large (in terms of
// <xsize, ysize, zsize>), photons scatter very often
// before exiting, so the code becomes incredibly slow
// without limiting the allowed number of scatterings
// (but doing THAT means photons will not explore the
// grid fully)
Grid::Grid(cuint _xres, cuint _yres, cuint _zres, cfloat _xsize, cfloat _ysize, cfloat _zsize, cfloat _den, cfloat _opa):
	xres(_xres), yres(_yres), zres(_zres),

	hxsize(_xsize * 0.5f), hysize(_ysize * 0.5f), hzsize(_zsize * 0.5f),
	 xsize(hxsize * 2.0f),  ysize(hysize * 2.0f),  zsize(hzsize * 2.0f),

	// vec3's for convenience
	hsize(hxsize, hysize, hzsize),
	fsize( xsize,  ysize,  zsize),

	density(_den),
	opacity(_opa) {

	Init();
}

void Grid::Init() {
	#ifdef VERBOSE
	std::cout << "[Grid::" << _FUNCTION__ << "]" << std::endl;
	std::cout << "\t(global) density: " << density << std::endl;
	std::cout << "\t(global) opacity: " << opacity << std::endl;
	std::cout << "\txres: " << xres << std::endl;
	std::cout << "\tyres: " << yres << std::endl;
	std::cout << "\tzres: " << zres << std::endl;
	std::cout << "\th-xsize: " << hsize.x << ", x-range: [" << -hsize.x << ", " << hsize.x << "]" << std::endl;
	std::cout << "\th-ysize: " << hsize.y << ", y-range: [" << -hsize.y << ", " << hsize.y << "]" << std::endl;
	std::cout << "\th-zsize: " << hsize.z << ", z-range: [" << -hsize.z << ", " << hsize.z << "]" << std::endl;
	#endif

	// number of faces in each dimension
	// is one larger than number of cells
	xfaces.resize(xres + 1, 0.0f);
	yfaces.resize(yres + 1, 0.0f);
	zfaces.resize(zres + 1, 0.0f);

	for (uint i = 0; i <= xres; i++) {
		xfaces[i] = (i * xsize) / xres;
	}
	for (uint i = 0; i <= yres; i++) {
		yfaces[i] = (i * ysize) / yres;
	}
	for (uint i = 0; i <= zres; i++) {
		zfaces[i] = (i * zsize) / zres;
	}


	// create the 3D density array (one element per cell)
	rhokappa.resize(xres, std::vector<fvec >(yres, fvec(zres, 0.0f)));

	// initialize and save the grid density array
	// std::fstream f0("GridDensity.dat", std::ios::out);
	// std::fstream f1("GridVolume.dat", std::ios::out);

	for (uint i = 0; i < xres; i++) {
		for (uint j = 0; j < yres; j++) {
			for (uint k = 0; k < zres; k++) {
				// constant-density cube
				if (i > uint(0.1f * xres) && i < uint(0.9f * xres) &&
					j > uint(0.1f * yres) && j < uint(0.9f * yres) &&
					k > uint(0.1f * zres) && k < uint(0.9f * zres)) {

					rhokappa[i][j][k] = density * opacity;
				} else {
					rhokappa[i][j][k] = 1.0f;
				}
			}

			// f0 << i << "\t" << j << "\t" <<  rk << endl;
		}
	}

	/*
	for (uint i = 0; i <= xres; i++) {
		for (uint j = 0; j <= yres; j++) {
			for (uint k = 0; k <= zres; k++) {
				f1 << xfaces[i] << "\t" << yfaces[j] << "\t" << zfaces[k] << std::endl;
			}
		}
	}
	*/

	// f0.close();
	// f1.close();
}



float Grid::GetMaxXDistance(const Photon& pho, const vec3& pcur) const {
	if (pho.dir.x > 0.0f) {
		return ((xsize - pcur.x) / pho.dir.x);
	} else if (pho.dir.x < 0.0f) {
		return (-pcur.x / pho.dir.x);
	} else {
		return (hsize.x * 10000.0f);
	}
}

float Grid::GetMaxYDistance(const Photon& pho, const vec3& pcur) const {
	if (pho.dir.y > 0.0f) {
		return ((ysize - pcur.y) / pho.dir.y);
	} else if (pho.dir.y < 0.0f) {
		return (-pcur.y / pho.dir.y);
	} else {
		return (hsize.y * 10000.0f);
	}
}

float Grid::GetMaxZDistance(const Photon& pho, const vec3& pcur) const {
	if (pho.dir.z > 0.0f) {
		return ((zsize - pcur.z) / pho.dir.z);
	} else if (pho.dir.z < 0.0f) {
		return (-pcur.z / pho.dir.z);
	} else {
		return (hsize.z * 10000.0f);
	}
}



// "taufind" (used to find the optical depth
// along a given photon direction for forced
// first scattering ONLY)
float Grid::GetMaxOpticalDepth(const Photon& p) const {
	float dpath = 0.0f;			// cumulative distance along photon path
	float dcell = 0.0f;			// distance to next cell wall (min(dx, dy, dz))
	float taurun = 0.0f;		// optical depth along photon path
	float taucell = 0.0f;

	vec3 pcur = p.pos + hsize;

	int celli = p.xcell;
	int cellj = p.ycell;
	int cellk = p.zcell;

	// maximum distance we can travel through grid in
	// x, y, z direction, distance to next face along
	// direction of travel
	cfloat dsx = GetMaxXDistance(p, pcur); float dx = 0.0f;
	cfloat dsy = GetMaxYDistance(p, pcur); float dy = 0.0f;
	cfloat dsz = GetMaxZDistance(p, pcur); float dz = 0.0f;

	// the largest distance we can travel through the grid is
	// the minimum of the x, y, and z-distances to each edge
	cfloat dmax = std::min(std::min(dsx, dsy), dsz);

	if (dmax < 0.01f) {
		// about to exit the grid
		return 0.0f;
	}


	while (dpath < dmax * 0.999f) {
		// get the optical depth to the next cell wall (ie.
		// "distance to wall of current cell * its opacity")
		taucell = dcell * rhokappa[celli][cellj][cellk];
		// if taurun + taucell < tau, then photon moves a
		// distance dcell (i.e. ends up on next cell wall)
		taurun += taucell;
		pcur += (p.dir * dcell);

		celli = GetXCellIdx(pcur.x - hsize.x);
		cellj = GetYCellIdx(pcur.y - hsize.y);
		cellk = GetZCellIdx(pcur.z - hsize.z);


		#ifdef DEBUG
		assert(IndicesInBounds(celli, cellj, cellk));
		#endif


		// find the distances dx, dy, and dz to the
		// x, y, and z cell-faces along the photon's
		// direction unit vector
		if (p.dir.x > 0.0f) {
			dx = (xfaces[celli + 1] - pcur.x) / p.dir.x;

			if (dx < EPSILON) {
				pcur.x = xfaces[celli + 1];
				celli += 1;
				dx = (xfaces[celli + 1] - pcur.x) / p.dir.x;
			}
		} else if (p.dir.x < 0.0f) {
			dx = (xfaces[celli] - pcur.x) / p.dir.x;

			if (dx < EPSILON) {
				pcur.x = xfaces[celli];
				dx = (xfaces[celli - 1] - pcur.x) / p.dir.x;
				celli -= 1;
			}
		} else if (p.dir.x == 0.0f) {
			dx = hsize.x * 10000.0f;
		}

		if (p.dir.y > 0.0f) {
			dy = (yfaces[cellj + 1] - pcur.y) / p.dir.y;

			if (dy < EPSILON) {
				pcur.y = yfaces[cellj + 1];
				cellj += 1;
				dy = (yfaces[cellj + 1] - pcur.y) / p.dir.y;
			}
		} else if (p.dir.y < 0.0f) {
			dy = (yfaces[cellj] - pcur.y) / p.dir.y;

			if (dy < EPSILON) {
				pcur.y = yfaces[cellj];
				dy = (yfaces[cellj - 1] - pcur.y) / p.dir.y;
				cellj -= 1;
			}
		} else if (p.dir.y == 0.0f) {
			dy = hsize.y * 10000.0f;
		}

		if (p.dir.z > 0.0f) {
			dz = (zfaces[cellk + 1] - pcur.z) / p.dir.z;

			if (dz < EPSILON) {
				pcur.z = zfaces[cellk + 1];
				cellk += 1;
				dz = (zfaces[cellk + 1] - pcur.z) / p.dir.z;
			}
		} else if (p.dir.z < 0.0f) {
			dz = (zfaces[cellk] - pcur.z) / p.dir.z;

			if (dz < EPSILON) {
				pcur.z = zfaces[cellk];
				dz = (zfaces[cellk - 1] - pcur.z) / p.dir.z;
				cellk -= 1;
			}
		} else if (p.dir.z == 0.0f) {
			dz = hsize.z * 10000.0f;
		}

		// distances are only zero if photon is on cell wall;
		// if so, then set them to arbitrary large values
		// since we will hit another wall
		if (dx == 0.0f || fabs(dx) < EPSILON) { dx = hsize.x * 10000.0f; }
		if (dy == 0.0f || fabs(dy) < EPSILON) { dy = hsize.y * 10000.0f; }
		if (dz == 0.0f || fabs(dz) < EPSILON) { dz = hsize.z * 10000.0f; }

		// find distance to next cell wall
		dcell = std::min(std::min(dx, dy), dz);

		if (dx < 0.0f) { dcell = std::min(dy, dz); }
		if (dy < 0.0f) { dcell = std::min(dx, dz); }
		if (dz < 0.0f) { dcell = std::min(dx, dy); }

		dpath += dcell;
	}

	return taurun;
}



void Grid::SetForcedFirstScatterPos(Photon& p, const RNG& rng) const {
	cfloat tEdgeMax  = GetMaxOpticalDepth(p);
	cfloat photonWgt = 1.0f - expf(-tEdgeMax);
	cfloat tRunMax   = -logf(1.0f - (rng.FUniformRandom() * photonWgt));

	Integrate(p, tRunMax, IM_FORCED);
}



// combines "tauint" (used to find the location of the
// forced first scattering event) and "tauint2" (used
// to find subsequent scattering positions within the
// grid from RANDOM optical depths)
bool Grid::Integrate(Photon& p, cfloat tau, IntegrateMode im) const {
	bool ret = true;

	float taurun = 0.0f, taucell = 0.0f;
	float dtmp = 0.0f, dpath = 0.0f, dcell = 0.0f;

	// translate by half the grid's represented dimensions
	// (so pcur corresponds to the photon cell indices)
	vec3 pcur = p.pos + hsize;

	int celli = p.xcell;
	int cellj = p.ycell;
	int cellk = p.zcell;


	#ifdef DEBUG
	assert(IndicesInBounds(celli, cellj, cellk));
	#endif


	cfloat dsx = GetMaxXDistance(p, pcur); float dx = 0.0f;
	cfloat dsy = GetMaxYDistance(p, pcur); float dy = 0.0f;
	cfloat dsz = GetMaxZDistance(p, pcur); float dz = 0.0f;

	cfloat dmax = std::min(std::min(dsx, dsy), dsz);

	if (dmax < .001f) {
		return false;
	}


	// start the path integration
	while (taurun < tau && dpath < dmax * 0.999f) {
		if (p.dir.x > 0.0f) {
			dx = (xfaces[celli + 1] - pcur.x) / p.dir.x;

			if (dx < EPSILON) {
				pcur.x = xfaces[celli + 1];
				celli += 1;
				dx = (xfaces[celli + 1] - pcur.x) / p.dir.x;
			}
		} else if (p.dir.x < 0.0f) {
			dx = (xfaces[celli] - pcur.x) / p.dir.x;

			if (dx < EPSILON) {
				pcur.x = xfaces[celli];
				dx = (xfaces[celli - 1] - pcur.x) / p.dir.x;
				celli -= 1;
			}
		} else if (p.dir.x == 0.0f) {
			dx = hsize.x * 10000.0f;
		}

		if (p.dir.y > 0.0f) {
			dy = (yfaces[cellj + 1] - pcur.y) / p.dir.y;

			if (dy < EPSILON) {
				pcur.y = yfaces[cellj + 1];
				cellj += 1;
				dy = (yfaces[cellj + 1] - pcur.y) / p.dir.y;
			}
		} else if (p.dir.y < 0.0f) {
			dy = (yfaces[cellj] - pcur.y) / p.dir.y;

			if (dy < EPSILON) {
				pcur.y = yfaces[cellj];
				dy = (yfaces[cellj - 1] - pcur.y) / p.dir.y;
				cellj -= 1;
			}
		} else if (p.dir.y == 0.0f) {
			dy = hsize.y * 10000.0f;
		}

		if (p.dir.z > 0.0f) {
			dz = (zfaces[cellk + 1] - pcur.z) / p.dir.z;

			if (dz < EPSILON) {
				pcur.z = zfaces[cellk + 1];
				cellk += 1;
				dz = (zfaces[cellk + 1] - pcur.z) / p.dir.z;
			}
		} else if (p.dir.z < 0.0f) {
			dz = (zfaces[cellk] - pcur.z) / p.dir.z;

			if (dz < EPSILON) {
				pcur.z = zfaces[cellk];
				dz = (zfaces[cellk - 1] - pcur.z) / p.dir.z;
				cellk -= 1;
			}
		} else if (p.dir.z == 0.0f) {
			dz = hsize.z * 10000.0f;
		}

		// distances can be zero if we
		// are right on a cell boundary
		if (dx == 0.0f || fabs(dx) < .001f) { dx = hsize.x * 10000.0f; }
		if (dy == 0.0f || fabs(dy) < .001f) { dy = hsize.y * 10000.0f; }
		if (dz == 0.0f || fabs(dz) < .001f) { dz = hsize.z * 10000.0f; }

		// find distance to next cell wall
		// (minimum value of dx, dy, and dz)
		dcell = std::min(std::min(dx, dy), dz);

		if (dx < 0.0f) { dcell = std::min(dy, dz); }
		if (dy < 0.0f) { dcell = std::min(dx, dz); }
		if (dz < 0.0f) { dcell = std::min(dx, dy); }

		// if moving through a "thicker" region of the grid
		// where the cells are denser, we reach tau quicker
		taucell = dcell * rhokappa[celli][cellj][cellk];

		// (if taurun + taucell > tau) then the photon's next
		// interaction event will be at distance d + dtmp, so
		// update photon position and cell indices; otherwise
		// photon moves over a distance dcell to the next cell
		// wall
		if (taurun + taucell >= tau) {
			dtmp = (tau - taurun) / rhokappa[celli][cellj][cellk];
			dpath += dtmp;
			taurun += taucell;

			pcur += (p.dir * dtmp);
		} else {
			dpath += dcell;
			taurun += taucell;

			pcur += (p.dir * dcell);
		}


		if (!GPosInBounds(pcur)) {
			// fight numerical inaccuracy
			ClampGPos(pcur);
		}


		celli = GetXCellIdx(pcur.x - hsize.x); // GetXCellIdx() already adds hsize.x to x!
		cellj = GetYCellIdx(pcur.y - hsize.y); // GetYCellIdx() already adds hsize.y to y!
		cellk = GetZCellIdx(pcur.z - hsize.z); // GetZCellIdx() already adds hsize.z to z!


		#ifdef DEBUG
		assert(IndicesInBounds(celli, cellj, cellk));
		#endif
	}


	switch (im) {
		case IM_RANDOM: {
			if (dpath >= dmax * 0.999f) {
				// photon escapes grid
				ret = false;
			} else {
				// photon stays inside, update
				// its pos and cell indices to
				// the next interaction coords
				p.UpdatePos(dpath);

				if (!RPosInBounds(p.pos)) {
					// fight numerical inaccuracy
					ClampRPos(p.pos);
				}

				p.UpdateCellIndices(*this);
			}
		} break;

		case IM_FORCED: {
			p.pos = pcur - hsize;

			if (!RPosInBounds(p.pos)) {
				// fight numerical inaccuracy
				ClampRPos(p.pos);
			}

			p.UpdateCellIndices(*this);
		} break;
	}

	return ret;
}



bool Grid::IndicesInBounds(cint i, cint j, cint k) const {
	cbool b0 = (i >=       0  && j >=       0  && k >=       0 );
	cbool b1 = (i < int(xres) && j < int(yres) && k < int(zres));

	return (b0 && b1);
}



// test a position in world-space [-hsize, hsize)
bool Grid::RPosInBounds(const vec3& p) const {
	cbool b0 = (p.x >= -hsize.x && p.x < hsize.x);
	cbool b1 = (p.y >= -hsize.y && p.y < hsize.y);
	cbool b2 = (p.z >= -hsize.z && p.z < hsize.z);
	return (b0 && b1 && b2);
}

// clamp a position in world-space
void Grid::ClampRPos(vec3& p) const {
	if (p.x < -hsize.x) { p.x = -hsize.x;           }
	if (p.x >= hsize.x) { p.x =  hsize.x - EPSILON; }

	if (p.y < -hsize.y) { p.y = -hsize.y;           }
	if (p.y >= hsize.y) { p.y =  hsize.y - EPSILON; }

	if (p.z < -hsize.z) { p.z = -hsize.z;           }
	if (p.z >= hsize.z) { p.z =  hsize.z - EPSILON; }
}


// test a position in grid-space [0, size)
bool Grid::GPosInBounds(const vec3& p) const {
	cbool b0 = (p.x >= 0.0f && p.x < xsize);
	cbool b1 = (p.y >= 0.0f && p.y < ysize);
	cbool b2 = (p.z >= 0.0f && p.z < zsize);
	return (b0 && b1 && b2);
}

// clamp a position in grid-space
void Grid::ClampGPos(vec3& p) const {
	if (p.x <   0.0f) { p.x =             0.0f; }
	if (p.x >= xsize) { p.x =  xsize - EPSILON; }

	if (p.y <   0.0f) { p.y =             0.0f; }
	if (p.y >= ysize) { p.y =  ysize - EPSILON; }

	if (p.z <   0.0f) { p.z =             0.0f; }
	if (p.z >= zsize) { p.z =  zsize - EPSILON; }
}
