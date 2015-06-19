#ifndef RNG_HDR
#define RNG_HDR

#include <cstdlib>
#ifdef THREADED
#include <boost/thread/mutex.hpp>
#else
#include <ctime>
#endif

#include "Constants.hpp"

struct RNG {
	RNG() {
		static bool inited = false;

		if (!inited) {
			inited = true;
			srandom(time(0x0));
		}
	}

	// returns a uniformly distributed
	// number in the range [0.0f, 1.0f]
	float FUniformRandom() const {
		#ifdef THREADED
		boost::mutex::scoped_lock l(rngMutex);
		#endif

		return (random() / FRAND_MAX);
	}

	#ifdef THREADED
	// needs to be mutable because FUR() is const
	mutable boost::mutex rngMutex;
	#endif
};

#endif
