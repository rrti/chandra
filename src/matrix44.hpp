#ifndef MATRIX44_HDR
#define MATRIX44_HDR

#include "vec3.hpp"

struct matrix44 {
	matrix44() {
		// load the identity matrix
		m[ 0] = m[ 5] = m[10] = m[15] = 1.0f;
		m[ 1] = m[ 2] = m[ 3] = m[ 4] = 0.0f;
		m[ 6] = m[ 7] = m[ 8] = m[ 9] = 0.0f;
		m[11] = m[12] = m[13] = m[14] = 0.0f;
	}

	inline void setCol(uint i, const vec4& v) {
		switch (i) {
			case 0: { m[ 0] = v.x; m[ 1] = v.y; m[ 2] = v.z; m[ 3] = v.w; } break;
			case 1: { m[ 4] = v.x; m[ 5] = v.y; m[ 6] = v.z; m[ 7] = v.w; } break;
			case 2: { m[ 8] = v.x; m[ 9] = v.y; m[10] = v.z; m[11] = v.w; } break;
			case 3: { m[12] = v.x; m[13] = v.y; m[14] = v.z; m[15] = v.w; } break;
		}
	}
	inline void setRow(uint i, const vec4& v) {
		switch (i) {
			case 0: { m[0] = v.x; m[4] = v.y; m[ 8] = v.z; m[12] = v.w; } break;
			case 1: { m[1] = v.x; m[5] = v.y; m[ 9] = v.z; m[13] = v.w; } break;
			case 2: { m[2] = v.x; m[6] = v.y; m[10] = v.z; m[14] = v.w; } break;
			case 3: { m[3] = v.x; m[7] = v.y; m[11] = v.z; m[15] = v.w; } break;
		}
	}

	inline vec4 mul(const vec4& v) const {
		//             | x
		//             | y
		//             | z
		//             | w
		// ------------+------------------
		// a  b  c  d  | ax + by + cz + dw
		// e  f  g  h  | ex + fy + gz + hw
		// i  j  k  l  | ix + jy + kz + lw
		// m  n  o  p  | mx + ny + oz + pw
		vec4 r;
		r.x = m[0] * v.x + m[4] * v.y + m[ 8] * v.z + m[12] * v.w;
		r.y = m[1] * v.x + m[5] * v.y + m[ 9] * v.z + m[13] * v.w;
		r.z = m[2] * v.x + m[6] * v.y + m[10] * v.z + m[14] * v.w;
		r.w = m[3] * v.x + m[7] * v.y + m[11] * v.z + m[15] * v.w;
		return r;
	}

	inline matrix44& mul(cfloat s) {
		m[ 0] *= s; m[ 1] *= s; m[ 2] *= s; m[ 3] *= s;
		m[ 4] *= s; m[ 5] *= s; m[ 6] *= s; m[ 7] *= s;
		m[ 8] *= s; m[ 9] *= s; m[10] *= s; m[11] *= s;
		m[12] *= s; m[13] *= s; m[14] *= s; m[15] *= s;
		return *this;
	}

	// elements are stored in
	// column-major order, ie.
	// m[0], ..., m[3] is the
	// first column
	float m[16];
};

#endif
