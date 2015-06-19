#ifndef MCRTSIM_HDR
#define MCRTSIM_HDR

#include <vector>

#ifdef THREADED
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#endif

// because of composition pattern, we
// need to include _everything_ here
// (or store pointers and dynamically
// allocate...)
#include "Constants.hpp"
#include "TypeDefs.hpp"
#include "ThreadData.hpp"
#include "SimData.hpp"
#include "RNG.hpp"
#include "Observer.hpp"
#include "Slab.hpp"
#include "Grid.hpp"
#include "Dust.hpp"
#include "Image.hpp"
#include "TDFStackParser.hpp"

struct MCRTSim {
	MCRTSim(CTDFStackParser& p, ProjectMode pm, ScatterMode sm, TraceMode tm);

	void WriteDebugData(const std::string& dir, const std::string& fname, const std::list<vec3>& data) const;
	void WriteGlobalStatsData(const std::string& dir);

	void InitGlobalStatsData();

	void MergeDataHelper(const ThreadData& td);

	#ifdef THREADED
	void MergeGlobalStatsData();
	#endif

	uint TracePhotonSlab(Photon& p, ThreadData* td, uint* pixelsHit);
	uint TracePhotonGrid(Photon& p, ThreadData* td, uint* pixelsHit);
	void TracePhotons(ThreadData* td);

	#ifdef THREADED
	void SpawnThreads(cuint numThreads);
	void JoinThreads();
	#endif

	void Run(cuint numThreads, cuint _numPhotons, cuint _numScatEvents, cbool verbose);

	void UnitTest0();
	void UnitTest1();
	void UnitTest2();



	RNG rng;
	Observer ob;

	Slab slab;
	Grid grid;
	Dust dust;

	#ifdef THREADED
	std::vector<boost::thread*> threads;
	std::vector<ThreadData> threadData;
	boost::mutex coutMutex;
	#endif

	SimData sd;
};

#endif
