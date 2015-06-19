#ifndef THREADDATA_HDR
#define THREADDATA_HDR

#include <list>
#include <vector>

#include "TypeDefs.hpp"
#include "vec3.hpp"
#include "Image.hpp"

// represents the "local data"
// of a photon-tracing thread
struct ThreadData {
	ThreadData(cuint tID, cuint tNumPhotons, cuint nBins, cuint nLvls, cuint imgW, cuint imgH):
		threadID(tID), numPhotons(tNumPhotons), img(imgW, imgH) {

		nrg.resize(nBins + 1, 0.0);
		err.resize(nBins + 1, 0.0);

		jPos.resize(nLvls, 0.0); jNeg.resize(nLvls, 0.0);
		hPos.resize(nLvls, 0.0); hNeg.resize(nLvls, 0.0);
		kPos.resize(nLvls, 0.0); kNeg.resize(nLvls, 0.0);
	}

	uint threadID;
	uint numPhotons;
	Image img;
	std::list<vec3> exitPosDistr;
	std::list<vec3> exitDirDistr;
	std::list<vec3> photonPosHst;
	fvec nrg;
	fvec err;
	fvec jPos, jNeg;
	fvec hPos, hNeg;
	fvec kPos, kNeg;
	std::list<uint> timings;
};

#endif
