#ifndef IMAGE_HDR
#define IMAGE_HDR

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Constants.hpp"
#include "TypeDefs.hpp"

struct Image {
	Image(cuint _w, cuint _h) {
		w = _w; hw = w >> 1;
		h = _h; hh = h >> 1;
		m =  0;

		d.resize(w, std::vector<uint>(h, 0));
	}

	void Write(const std::string& dir, const std::string& fname) const {
		std::cout << "[Image::" << __FUNCTION__ << "]" << std::endl;

		std::string path = dir + "/" + fname;
		std::ofstream f(path.c_str(), std::ios::out);
		std::stringstream ss;

		ss << "P3" << std::endl;
		ss << "## Plane-Parallel *tropic Scattering" << std::endl;
		ss << w << " " << h << std::endl;
		ss << m << std::endl;

		uint zps = 0;
		cuint s = w * h;
		cfloat maxEnergy = logf(m);

		// output in row-major order (all
		// columns in row 0, then all in
		// row 1, etc)
		for (uint i = 0; i < s; i++) {
			cuint ix = i % w;
			cuint iy = i / w;
			cuint px = d[ix][iy];

			// ss << px << " " << px << " " << px << endl;
			ss << uint((log(d[ix][iy]) / maxEnergy) * PXL_MAX_VAL) << " " << px << " " << px << std::endl;

			if (px == 0) {
				zps++;
			}
		}

		f << ss.str();
		f.close();

		#ifdef VERBOSE
		std::cout << "\tnumber of zero-pixels: " << zps;
		std::cout << " (total: " << s << ")" << std::endl;
		#endif
	}

	uint w, hw;		// width, half-width in pixels
	uint h, hh;		// height, half-height in pixels
	uint m;			// maximum pixel-value

	std::vector<std::vector<uint> > d;
};

#endif
