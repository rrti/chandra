#ifndef SIMDATA_HDR
#define SIMDATA_HDR

#include <iostream>
#include <list>

#include "vec3.hpp"
#include "TypeDefs.hpp"
#include "Image.hpp"

struct SimData {
	SimData(
		ProjectMode pm,
		ScatterMode sm,
		TraceMode tm,
		const std::string& ddir,
		cuint nPhotons,
		cuint nScatEvts,
		cuint nMuBins,
		cuint imgW,
		cuint imgH
	):
		numPhotons(nPhotons),
		numScatEvents(nScatEvts),
		numMuBins(nMuBins),
		dataDir(ddir),
		verbose(false),
		projectMode(pm),
		scatterMode(sm),
		traceMode(tm),
		img(imgW, imgH) {
	}

	uint numPhotons;
	uint numScatEvents;
	uint numMuBins;
	std::string dataDir;
	bool verbose;

	ProjectMode projectMode;		// perspective or orthographic
	ScatterMode scatterMode;		// iso- or anisotropic
	TraceMode traceMode;			// slab or grid

	Image img;						// the simulation's visual output (radiant energy)

	// intensity moments
	fvec jPos, jNeg;
	fvec hPos, hNeg;
	fvec kPos, kNeg;

	// per-bin total energy, intensity, deviation, and error
	fvec energy;
	fvec intens;
	fvec sigmaI;
	fvec errorI;
	fvec thetas;

	std::list<vec3> exitPosDistr;
	std::list<vec3> exitDirDistr;
	std::list<vec3> photonPosHst;

	std::list<uint> timings;
};

#endif
