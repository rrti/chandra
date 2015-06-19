#include <iostream>
#include <fstream>
#include <sstream>

#include <cmath>
#include <cstdlib>
#include <sys/stat.h>

#include <SDL/SDL_timer.h>

#include "MCRTSim.hpp"
#include "vec3.hpp"
#include "Photon.hpp"

using std::cout;
using std::endl;
using std::min;
using std::max;

MCRTSim::MCRTSim(CTDFStackParser& p, ProjectMode pm, ScatterMode sm, TraceMode tm):
	ob(
		p.GetVec("root\\observer", "pos"),
		p.GetValue<float>("root\\observer", "the"), p.GetValue<float>("root\\observer", "phi")
	),
	slab(
		p.GetValue<float>("root\\slab", "zmin"), p.GetValue<float>("root\\slab", "zmax"),
		p.GetValue<float>("root\\slab", "density"), p.GetValue<float>("root\\slab", "opacity"),
		p.GetValue<float>("root\\slab", "numlevels")
	),
	grid(
		p.GetValue<uint>( "root\\grid", "xres" ), p.GetValue<uint>( "root\\grid", "yres" ), p.GetValue<uint>( "root\\grid", "zres" ),
		p.GetValue<float>("root\\grid", "xsize"), p.GetValue<float>("root\\grid", "ysize"), p.GetValue<float>("root\\grid", "zsize"),
		p.GetValue<float>("root\\grid", "density"), p.GetValue<float>("root\\grid", "opacity")
	),
	dust(
		p.GetValue<float>("root\\dust", "a"),
		p.GetValue<float>("root\\dust", "g"),
		p.GetValue<float>("root\\dust", "pl"),
		p.GetValue<float>("root\\dust", "pc"),
		p.GetValue<float>("root\\dust", "sc")
	),
	sd(
		pm,
		sm,
		tm,
		p.GetValue<std::string>("root\\simulation", "datadir"),
		p.GetValue<uint>("root\\simulation", "numphotons"),
		p.GetValue<uint>("root\\simulation", "numscatevents"),
		p.GetValue<uint>("root\\simulation", "nummubins"),
		p.GetValue<uint>("root\\image", "width"), p.GetValue<uint>("root\\image", "height")
	) {
}

void MCRTSim::WriteDebugData(const std::string& dir, const std::string& fname, const std::list<vec3>& data) const {
	std::string path = dir + "/" + fname;
	std::ofstream f(path.c_str(), std::ios::out);
	std::stringstream ss;

	for (std::list<vec3>::const_iterator it = data.begin(); it != data.end(); it++) {
		ss << ((*it).x) << "\t" << ((*it).y) << "\t" << ((*it).z) << endl;
	}

	f << ss.str();
	f.close();
}

void MCRTSim::WriteGlobalStatsData(const std::string& dir) {
	mkdir(dir.c_str(), S_IRWXU | S_IRWXG);

	sd.img.Write(dir, "Image.ppm");

	char buf0[64] = {'\0'}; snprintf(buf0, 63, "%s/J-moment.dat",         dir.c_str());
	char buf1[64] = {'\0'}; snprintf(buf1, 63, "%s/H-moment.dat",         dir.c_str());
	char buf2[64] = {'\0'}; snprintf(buf2, 63, "%s/K-moment.dat",         dir.c_str());
	char buf3[64] = {'\0'}; snprintf(buf3, 63, "%s/ThetaIntensity.dat",   dir.c_str());
	char buf4[64] = {'\0'}; snprintf(buf4, 63, "%s/ThetaSigmaEnergy.dat", dir.c_str());
	char buf5[64] = {'\0'}; snprintf(buf5, 63, "%s/ExecutionTimings.dat", dir.c_str());

	char buf6[64] = {'\0'}; snprintf(buf6, 63, "%s/Moments.plt",          dir.c_str());
	char buf7[64] = {'\0'}; snprintf(buf7, 63, "%s/TheIns.plt",           dir.c_str());
	char buf8[64] = {'\0'}; snprintf(buf8, 63, "%s/TheSig.plt",           dir.c_str());
	char buf9[64] = {'\0'}; snprintf(buf9, 63, "%s/ExecutionTimings.plt", dir.c_str());

	std::fstream f0(buf0, std::ios::out);
	std::fstream f1(buf1, std::ios::out);
	std::fstream f2(buf2, std::ios::out);
	std::fstream f3(buf3, std::ios::out);
	std::fstream f4(buf4, std::ios::out);
	std::fstream f5(buf5, std::ios::out);

	std::fstream f6(buf6, std::ios::out);
	std::fstream f7(buf7, std::ios::out);
	std::fstream f8(buf8, std::ios::out);
	std::fstream f9(buf9, std::ios::out);


	for (uint i = 0; i < slab.numLevels; i++) {
		cfloat tau = slab.tauMax - ((slab.tauMax / slab.numLevels) * i);

		sd.jPos[i] /= sd.numPhotons;
		sd.jNeg[i] /= sd.numPhotons;
		sd.hPos[i] /= sd.numPhotons;
		sd.hNeg[i] /= sd.numPhotons;
		sd.kPos[i] /= sd.numPhotons;
		sd.kNeg[i] /= sd.numPhotons;

		f0 << tau << "\t" << (sd.jPos[i] - 0.0 * fabs(sd.jNeg[i])) << endl;
		f1 << tau << "\t" << (sd.hPos[i] - 1.0 * fabs(sd.hNeg[i])) << endl; 
		f2 << tau << "\t" << (sd.kPos[i] - 0.0 * fabs(sd.kNeg[i])) << endl;
	}

	f0.close();
	f1.close();
	f2.close();


	// float totalNrg = 0.0;
	// float totalInt = 0.0;

	for (uint i = 0; i < sd.numMuBins; i++) {
		if (sd.energy[i] > 0.0) {
			cfloat nrgNorm = sd.numPhotons * 2 * cosf(sd.thetas[i] * DEG2RAD);

			// energy[k] and errorI[k] are increased by 1 when
			// a photon exits the slab on the observer's side;
			// all the other measures are derived
			//
			//  E_i = N_i / N_t
			// sE_i = E_i / sqrt(N_i)
			// sE_i = N_i / N_t / sqrt(N_i)
			// sE_i = N_i / (N_t * sqrt(N_i))
			sd.intens[i] = sd.energy[i] / nrgNorm * sd.numMuBins;
			sd.sigmaI[i] = (sd.energy[i] / sd.numPhotons) / sqrtf(sd.energy[i]);
			sd.energy[i] = sd.energy[i] / sd.numPhotons;

			// totalNrg += sd.energy[i];
			// totalInt += sd.intens[i];
		}

		f3 << sd.thetas[i] << "\t" << sd.intens[i] << endl;
		f4 << sd.thetas[i] << "\t" << sd.sigmaI[i] << endl;
	}

	f3.close();
	f4.close();


	cuint timingInt = sd.numPhotons / 10;
	uint timingPoint = timingInt;

	for (std::list<uint>::const_iterator git = sd.timings.begin(); git != sd.timings.end(); git++) {
		f5 << timingPoint << "\t" << *git << endl;
		timingPoint += timingInt;
	}

	f5.close();


//	f6 << "set term png" << endl;
//	f6 << "set output IntensityMoments.png" << endl;
	f6 << "set term x11" << endl;
	f6 << "set xlabel \"tau\"" << endl;
	f6 << "set ylabel \"IntensityMoments(tau)\"" << endl;
	f6 << "plot \'J-moment.dat\' with lines title \'J\',\\" << endl;
	f6 << "     \'H-moment.dat\' with lines title \'H\',\\" << endl;
	f6 << "     \'K-moment.dat\' with lines title \'K\'   " << endl;
	f6 << "pause -1" << endl;
	f6.close();

	f7 << "set term x11" << endl;
	f7 << "set xlabel \"theta\"" << endl;
	f7 << "set ylabel \"Intensity(theta)\"" << endl;
	f7 << "plot \'ThetaIntensity.dat\' with points pointtype 5" << endl;
	f7 << "pause -1" << endl;
	f7.close();

	f8 << "set term x11" << endl;
	f8 << "set xlabel \"theta\"" << endl;
	f8 << "set ylabel \"SigmaEnergy(theta)\"" << endl;
	f8 << "plot \'ThetaSigmaEnergy.dat\' with points pointtype 5" << endl;
	f8 << "pause -1" << endl;
	f8.close();

	f9 << "set term x11" << endl;
	f9 << "set xlabel \"Number of Photons\"" << endl;
	f9 << "set ylabel \"Execution Time (ms)\"" << endl;
	f9 << "plot \'ExecutionTimings.dat\' with lines" << endl;
	f9 << "pause -1" << endl;
	f9.close();
}

void MCRTSim::InitGlobalStatsData() {
	sd.jPos.resize(slab.numLevels, 0.0); sd.jNeg.resize(slab.numLevels, 0.0);
	sd.hPos.resize(slab.numLevels, 0.0); sd.hNeg.resize(slab.numLevels, 0.0);
	sd.kPos.resize(slab.numLevels, 0.0); sd.kNeg.resize(slab.numLevels, 0.0);

	// cos(theta) can be 1, so one
	// extra bin index is needed
	sd.energy.resize(sd.numMuBins + 1, 0.0); // written in TracePhotonSlab() only
	sd.intens.resize(sd.numMuBins,     0.0); // not read or written until WriteGlobalStatsData()
	sd.sigmaI.resize(sd.numMuBins,     0.0); // not read or written until WriteGlobalStatsData()
	sd.errorI.resize(sd.numMuBins + 1, 0.0); // written in TracePhotonSlab() only
	sd.thetas.resize(sd.numMuBins,     0.0); // becomes read-only after initialization

	// per-bin delta-theta and half-width
	cfloat dt = 1.0 / sd.numMuBins;
	cfloat hw = dt * 0.5;

	for (uint i = 1; i <= sd.numMuBins; i++) {
		sd.thetas[i - 1] = acosf((i - 1) * dt + hw) * RAD2DEG;
	}
}


void MCRTSim::MergeDataHelper(const ThreadData& td) {
	std::list<vec3>::const_iterator lit;

	for (lit = td.exitPosDistr.begin(); lit != td.exitPosDistr.end(); lit++) {
		sd.exitPosDistr.push_back(*lit);
	}
	for (lit = td.exitDirDistr.begin(); lit != td.exitDirDistr.end(); lit++) {
		sd.exitDirDistr.push_back(*lit);
	}
	for (lit = td.photonPosHst.begin(); lit != td.photonPosHst.end(); lit++) {
		sd.photonPosHst.push_back(*lit);
	}

	for (uint w = 0; w < sd.img.w; w++) {
		for (uint h = 0; h < sd.img.h; h++) {
			sd.img.d[w][h] += td.img.d[w][h];
			sd.img.d[w][h] = min(PXL_MAX_VAL, sd.img.d[w][h]);
			sd.img.m = max(sd.img.m, sd.img.d[w][h]);
		}
	}

	for (uint j = 0; j < sd.numMuBins; j++) {
		sd.energy[j] += td.nrg[j];
		sd.errorI[j] += td.err[j];
	}

	for (uint j = 0; j < slab.numLevels; j++) {
		sd.jPos[j] += td.jPos[j]; sd.jNeg[j] += td.jNeg[j];
		sd.hPos[j] += td.hPos[j]; sd.hNeg[j] += td.hNeg[j];
		sd.kPos[j] += td.kPos[j]; sd.kNeg[j] += td.kNeg[j];
	}


	if (sd.timings.empty()) {
		std::list<uint>::const_iterator tit;

		for (tit = td.timings.begin(); tit != td.timings.end(); tit++) {
			sd.timings.push_back(*tit);
		}
	} else {
		std::list<uint>::iterator git = sd.timings.begin();
		std::list<uint>::const_iterator tit;

		for (tit = td.timings.begin(); tit != td.timings.end(); tit++) {
			// single-core CPU:
			//   1 thread:
			//       total program running time is T
			//   4 threads:
			//       total program running time is T
			//       per-thread running time is either
			//       T/4 (AAAA BBBB CCCC DDDD schedule)
			//       or T (ABCD ABCD ABCD ABCD schedule)
			// quad-core CPU:
			//   1 thread:
			//       total program running time is T
			//   4 threads:
			//       total program running time is T/4
			//       per-thread running time is T/4
			//
			// with multiple threads, the total running time
			// at each measure point is _not_ the sum of the
			// per-thread running times!
			*git = max(*git, *tit); git++;
		}
	}
}

#ifdef THREADED
void MCRTSim::MergeGlobalStatsData() {
	// merge the per-thread "raw" data
	for (uint tID = 0; tID < threads.size(); tID++) {
		MergeDataHelper(threadData[tID]);
	}
}
#endif



uint MCRTSim::TracePhotonSlab(Photon& p, ThreadData* td, uint* pixelsHit) {
//	td->photonPosHst.push_back(p.pos);

	while (slab.PosInBounds(p.pos)) {
		const vec3 opos = p.pos;

		p.UpdatePos(slab.GetRandomL(rng));
//		td->photonPosHst.push_back(p.pos);

		slab.UpdateIntensityMoments(td, opos, p);

		if (!slab.PosInBounds(p.pos)) {
//			if (p.pos.z > slab.zmax) {
//				td->exitPosDistr.push_back(p.pos);
//				td->exitDirDistr.push_back(p.dir);
//			}

			break;
		}

		if (rng.FUniformRandom() < dust.a) {
			p.Scatter(rng, dust, sd.scatterMode);
		} else {
			// photon absorbed into "nebula"
			p.absorbed = true;
			break;
		}
	}

	if (!p.absorbed) {
		if (p.pos.z < slab.zmin) {
			// photon exited on wrong side,
			// force emission of a new one
			return 0;
		} else {
			switch (sd.projectMode) {
				case PM_PERSPECTIVE: { (*pixelsHit) += p.PerspProject(ob, &(td->img)); } break;
				case PM_ORTHOGRAPHIC: { (*pixelsHit) += p.OrthoProject(ob, &(td->img)); } break;
			}

			cuint bin = uint(sd.numMuBins * fabs(p.cthe));
			td->nrg[bin] += 1.0;
			td->err[bin] += 1.0;

			return 1;
		}
	}

	// photon was absorbed, continue
	// with the next one (no respawn)
	return 1;
}

uint MCRTSim::TracePhotonGrid(Photon& p, ThreadData* td, uint* pixelsHit) {
	grid.SetForcedFirstScatterPos(p, rng);
	p.Scatter(rng, dust, sd.scatterMode);

	// NOTE: if we limit the number of scatterings then we
	// should also implement weight "peeling" (of direct-
	// and scatter-photons)?
	while (grid.Integrate(p, -logf(rng.FUniformRandom() + 0.01f), IM_RANDOM)) {
		if (rng.FUniformRandom() < dust.a) {
			p.Scatter(rng, dust, sd.scatterMode);

			if (p.numScatEvents >= sd.numScatEvents) {
				break;
			}
		} else {
			p.absorbed = true;
			break;
		}
	}

	if (!p.absorbed) {
		if (ob.zdir.dot3D(p.dir) >= 0.0f) {
			// photon direction points away from the
			// image plane, so can never hit a pixel
			// (NOTE: technically this doesn't mean
			// we shouldn't increase the counter, but
			// we return 0 here for consistency with
			// TracePhotonSlab())
			return 0;
		} else {
			switch (sd.projectMode) {
				case PM_PERSPECTIVE: { (*pixelsHit) += p.PerspProject(ob, &(td->img)); } break;
				case PM_ORTHOGRAPHIC: { (*pixelsHit) += p.OrthoProject(ob, &(td->img)); } break;
			}

			return 1;
		}
	}

	return 1;
}

// independently executed by each thread
void MCRTSim::TracePhotons(ThreadData* td) {
	uint photonNum = 1;
	uint pixelsHit = 0;
	bool newPhoton = false;

	cuint ticks = SDL_GetTicks();
	cuint timingInt = td->numPhotons / 10;
	cuint msgInt = 1000;

	// TODO: make this configurable
	// (and allow for more than one)
	const vec3 sourcePos(0.0f, 0.0f, 0.0f);

	// calculate the scattered-light contribution
	while (photonNum <= td->numPhotons) {
		Photon p(rng, sourcePos, grid);

		switch (sd.traceMode) {
			case TM_SLAB: { newPhoton = (TracePhotonSlab(p, td, &pixelsHit) > 0); } break;
			case TM_GRID: { newPhoton = (TracePhotonGrid(p, td, &pixelsHit) > 0); } break;
		}

		if (newPhoton) {
			photonNum++;

			if ((photonNum % timingInt) == 0) {
				td->timings.push_back(SDL_GetTicks() - ticks);
			}

			if (sd.verbose && (photonNum % msgInt) == 0) {
				#ifdef THREADED
				boost::mutex::scoped_lock l(coutMutex);
				#endif

				cout << "[MCRTSim::" << __FUNCTION__ << "] [thread " << td->threadID << "] ";
				cout << photonNum << " photons traced";
				cout << " (" << (td->numPhotons - photonNum);
				cout << " remaining)" << endl;
			}
		}
	}

	if (sd.verbose) {
		#ifdef THREADED
		boost::mutex::scoped_lock l(coutMutex);
		#endif

		cout << "[MCRTSim::" << __FUNCTION__ << "] [thread " << td->threadID << "]";
		cout << " number of photons that struck an image pixel: " << pixelsHit;
		cout << endl;
	}
}



#ifdef THREADED
void MCRTSim::SpawnThreads(cuint numThreads) {
	threads.resize(numThreads, 0x0);

	if (sd.numPhotons != ((sd.numPhotons / numThreads) * numThreads)) {
		// for now, numThreads must evenly divide numPhotons
		return;
	}

	for (uint tID = 0; tID < numThreads; tID++) {
		threadData.push_back(ThreadData(tID, (sd.numPhotons / numThreads), sd.numMuBins, slab.numLevels, sd.img.w, sd.img.h));
	}

	for (uint tID = 0; tID < numThreads; tID++) {
		threads[tID] = new boost::thread(boost::bind(&MCRTSim::TracePhotons, this, &threadData[tID]));
	}
}

void MCRTSim::JoinThreads() {
	for (uint tID = 0; tID < threads.size(); tID++) {
		threads[tID]->join();
		delete threads[tID];
		threads[tID] = 0x0;
	}
}
#endif

void MCRTSim::Run(cuint numThreads, cuint _numPhotons, cuint _numScatEvents, cbool _verbose) {
	InitGlobalStatsData();

	sd.numPhotons = _numPhotons;
	sd.numScatEvents = _numScatEvents;
	sd.verbose = _verbose;

	#ifdef THREADED
	SpawnThreads(numThreads);
	JoinThreads();
	MergeGlobalStatsData();
	#else
	ThreadData td(0, sd.numPhotons, sd.numMuBins, slab.numLevels, sd.img.w, sd.img.h);
	TracePhotons(&td);
	MergeDataHelper(td);
	#endif

	WriteGlobalStatsData(sd.dataDir);
}






void MCRTSim::UnitTest0() {
	cfloat g  = dust.g;
	cfloat gg = dust.gg;
	cfloat pl = dust.pl;
	cfloat pc = dust.pc;
	cfloat sc = dust.sc;
	cfloat xi = 0.01f;

	cfloat n = (1.0f - gg) / (1.0f - g + 2.0f * g * xi);
	cfloat cTHE = max(-1.0f, min(1.0f, (1.0f / (2.0f * g))   * ((1.0f + gg) - (n * n))));
	cfloat cTHEsq = cTHE * cTHE;

	cfloat  THE  = acosf(cTHE);
	cfloat  THEf = THE * (1.0f + 3.13f * sc * expf((-7.0f * THE) / M_PI));
	cfloat cTHEf = cosf(THEf);
	// cfloat  PHI  = M_PI * (2.0f * rng.FUniformRandom() - 1.0f);
	// cfloat cPHI  = cosf(PHI);

	float p1 = 0.0f;
	float p2 = 0.0f;
	float p3 = 0.0f;
	float p4 = 0.0f;

	dust.GetMatrixParameters(p1, p2, p3, p4, cTHE, cTHEsq);

	cfloat P1 = ((1.0f - gg) / powf((1.0f + gg - 2.0f * g * cTHE), 1.5f));
	cfloat P2 = (-pl * P1 * ((1.0f - cTHE * cTHE) / (1.0f + cTHE * cTHE)));
	cfloat P3 = (P1 * ((2.0f * cTHE) / (1.0f + cTHE * cTHE)));
	cfloat P4 = (-pc * P1 * ((1.0f - cTHEf * cTHEf) / (1.0f + cTHEf * cTHEf)));

	cout << "[MCRTSim::" << __FUNCTION__ << "]" << endl;
	cout << "\tTHE: " << THE << ", cTHE: " << cTHE << ", THEf: " << THEf << endl;
	cout << "\tp1: " << p1 << ", p2: " << p2 << ", p3: " << p3 << ", p4: " << p4 << endl;
	cout << "\tP1: " << P1 << ", P2: " << P2 << ", P3: " << P3 << ", P4: " << P4 << endl;

	Photon p(rng, vec3(0.0f, 0.0f, 0.0f), grid);

	p.UpdateAngles(  0.0f * DEG2RAD, 0.0f * DEG2RAD); p.UpdateDir(); std::cout << p.str() << std::endl;
	p.UpdateAngles( 90.0f * DEG2RAD, 0.0f * DEG2RAD); p.UpdateDir(); std::cout << p.str() << std::endl;
	p.UpdateAngles(180.0f * DEG2RAD, 0.0f * DEG2RAD); p.UpdateDir(); std::cout << p.str() << std::endl;

	printf("\n");

	p.UpdateAngles(  0.0f * DEG2RAD, 90.0f * DEG2RAD); p.UpdateDir(); std::cout << p.str() << std::endl;
	p.UpdateAngles( 90.0f * DEG2RAD, 90.0f * DEG2RAD); p.UpdateDir(); std::cout << p.str() << std::endl;
	p.UpdateAngles(180.0f * DEG2RAD, 90.0f * DEG2RAD); p.UpdateDir(); std::cout << p.str() << std::endl;

	printf("\n");

	p.UpdateAngles(0.0f * DEG2RAD,   0.0f * DEG2RAD); p.UpdateDir(); std::cout << p.str() << std::endl;
	p.UpdateAngles(0.0f * DEG2RAD,  90.0f * DEG2RAD); p.UpdateDir(); std::cout << p.str() << std::endl;
	p.UpdateAngles(0.0f * DEG2RAD, 180.0f * DEG2RAD); p.UpdateDir(); std::cout << p.str() << std::endl;
}

void MCRTSim::UnitTest1() {
	cout << "[MCRTSim::" << __FUNCTION__ << "]" << endl;

	// Photon p0(rng, vec3(22.73f, 3.47f, 50.00f - EPSILON), grid);

	while (true) {
		Photon p(rng, vec3(0.0f, 0.0f, 0.0f), grid);
		p.dir = vec3(0.0f, 0.0f, 1.0f);

		cfloat dsx = grid.GetMaxXDistance(p, p.pos + grid.hsize);
		cfloat dsy = grid.GetMaxYDistance(p, p.pos + grid.hsize);
		cfloat dsz = grid.GetMaxZDistance(p, p.pos + grid.hsize);

		cout << "\tdsx: " << dsx << ", x-cell idx:" << grid.GetXCellIdx(p.pos.x) << endl;
		cout << "\tdsy: " << dsy << ", y-cell idx:" << grid.GetYCellIdx(p.pos.y) << endl;
		cout << "\tdsz: " << dsz << ", z-cell idx:" << grid.GetZCellIdx(p.pos.z) << endl;

		cout << "\tspawned grid photon, pos: " << p.pos.str() << ", dir: " << p.dir.str() << ", max. opt. depth: " << grid.GetMaxOpticalDepth(p) << endl;

		while (grid.Integrate(p, -logf(rng.FUniformRandom()), IM_RANDOM)) {
			cout << "\t\tp.pos: " << p.pos.str() << endl;
			cout << "\t\tp.dir: " << p.dir.str() << endl;
			cout << "\t\tp.pol: " << p.pol.x << ", " << p.pol.y << ", " << p.pol.z << ", " << p.pol.w << endl;
			cout << "\t\tp.xcell: " << p.xcell << ", p.ycell: " << p.ycell << ", p.zcell: " << p.zcell << endl;
			cout << endl;

			// no scattering here so we can easily track when photon leaves grid
			// p.Scatter(rng, dust, scatterMode);
		}
	}

	/*
	// cTHE should be isotropically distributed for g close
	// to 0 and very anisotropically for g close to 1 ==>
	// confirmed, scattering assymetry is reproduced
	std::fstream f("AnisoScatDirDist.dat", std::ios::out);

	for (uint i = 0; i < 12345U; i++) {
		cfloat g  = dust.g;
		cfloat gg = dust.gg;
		cfloat xi = rng.FUniformRandom();

		cfloat n = (1.0f - gg) / (1.0f - g + 2.0f * g * xi);
		cfloat cTHE = min(1.0f, (1.0f / (2.0f * g))   * ((1.0f + gg) - (n * n)));

		f << cTHE << endl;
	}

	f.close();
	*/
}

