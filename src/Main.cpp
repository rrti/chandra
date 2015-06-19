#include "MCRTSim.hpp"

int main(int argc, char** argv) {
	CTDFStackParser p("params.tdf");

	uint numThreads         =  p.GetValue<uint>("root\\simulation", "numthreads");
	uint numPhotons         =  p.GetValue<uint>("root\\simulation", "numphotons");
	uint numScatEvents      =  p.GetValue<uint>("root\\simulation", "numscatevents");

	ProjectMode projectMode = (p.GetValue<uint>("root\\simulation", "projectmode") == PM_PERSPECTIVE)? PM_PERSPECTIVE: PM_ORTHOGRAPHIC;
	ScatterMode scatterMode = (p.GetValue<uint>("root\\simulation", "scattermode") == SM_ISOTROPIC  )? SM_ISOTROPIC:   SM_ANISOTROPIC;
	TraceMode traceMode     = (p.GetValue<uint>("root\\simulation", "tracemode"  ) == TM_SLAB       )? TM_SLAB:        TM_GRID;

	bool enableVerbosity    = p.GetValue<bool>("root\\simulation", "enableverbosity");
	bool enableUnitTests    = p.GetValue<bool>("root\\simulation", "enableunittests");

	if (argc > 1) {
		for (int a = 0; a < argc; a++) {
			const std::string s(argv[a]);
			const bool moreArgs = (a < argc - 1);

			if (s == "--numthreads" && moreArgs) {
				numThreads = atoi(argv[a + 1]);
				continue;
			}

			if (s == "--numphotons" && moreArgs) {
				numPhotons = atoi(argv[a + 1]);
				continue;
			}

			if (s == "--numscatterevents" && moreArgs) {
				numScatEvents = atoi(argv[a + 1]);
				continue;
			}

			if (s == "--orthographic") {
				projectMode = PM_ORTHOGRAPHIC;
				continue;
			}

			if (s == "--anisotropic") {
				scatterMode = SM_ANISOTROPIC;
				continue;
			}

			if (s == "--grid") {
				traceMode = TM_GRID;
				continue;
			}

			if (s == "--test") {
				enableUnitTests = true;
				continue;
			}

			if (s == "--verbose") {
				enableVerbosity = true;
				continue;
			}
		}
	}

	if (enableVerbosity) {
		std::cout << "[" << __FUNCTION__ << "]" << std::endl;
		std::cout << "\tusing " << numThreads << " photon-trace " << ((numThreads > 1)? "threads": "thread");
		std::cout << std::endl;
		std::cout << "\tusing " << ((scatterMode == SM_ISOTROPIC)? "isotropic": "anisotropic");
		std::cout << " photon scattering mode";
		std::cout << std::endl;
		std::cout << "\tusing " << ((traceMode == TM_SLAB)? "slab-type": "grid-type");
		std::cout << " scattering atmosphere";
		std::cout << std::endl;
		std::cout << "\tusing " << ((projectMode == PM_PERSPECTIVE)? "perspective": "orthographic");
		std::cout << " image projection mode";
		std::cout << std::endl;
		std::cout << "\tnumber of photons to track: " << numPhotons << std::endl;
		std::cout << "\tmax. scatter-events per photon (while in grid): "  << numScatEvents << std::endl;
	}

	MCRTSim sim(p, projectMode, scatterMode, traceMode);

	if (enableUnitTests) {
		sim.UnitTest0();
		sim.UnitTest1();
	} else {
		sim.Run(numThreads, numPhotons, numScatEvents, enableVerbosity);
	}

	return 0;
}

