#ifndef OBSERVER_HDR
#define OBSERVER_HDR

#include <iostream>
#include "vec3.hpp"

struct Observer {
	Observer(const vec3& _pos, cfloat /*_the*/, cfloat /*_phi*/): pos(_pos) /*, the(_the), phi(_phi)*/ {
		// cthe = cosf(the); sthe = sinf(the);
		// cphi = cosf(phi); sphi = sinf(phi);

		// note: assumes the observer always
		// looks at the origin of the system
		dst = pos.len3D();

		zdir = -(pos.norm());
		xdir = zdir.cross(vec3(0.0f, 1.0f, 0.0f)).norm();
		ydir = -(zdir.cross(xdir).norm());

		#ifdef VERBOSE
		std::cout << "[" << __FUNCTION__ << "]" << std::endl;
		std::cout << "\tpos:   " << pos.str() << std::endl;
		std::cout << "\tdist:  " << dst << std::endl;
		std::cout << "\tx-dir: " << xdir.str() << std::endl;
		std::cout << "\ty-dir: " << ydir.str() << std::endl;
		std::cout << "\tz-dir: " << zdir.str() << std::endl;
		std::cout << std::endl;
		#endif
	}

	vec3 pos;
	vec3 xdir;
	vec3 ydir;
	vec3 zdir;

	float dst;

	// float the, phi;
	// float cthe, sthe;
	// float cphi, sphi;
};

#endif
