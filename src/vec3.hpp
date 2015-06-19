#ifndef VEC3_HDR
#define VEC3_HDR

#include <cmath>
#include <cstdio>
#include <iostream>

#include "Constants.hpp"
#include "TypeDefs.hpp"

struct vec3 {
	vec3(float _x = 0.0f, float _y = 0.0f, float _z = 0.0f): x(_x), y(_y), z(_z) {}
	vec3(const vec3& v): x(v.x), y(v.y), z(v.z) {}

	inline bool operator == (const vec3& v) {
		return (fabsf(x - v.x) < EPSILON && fabsf(y - v.y) < EPSILON && fabsf(z - v.z) < EPSILON);
	}

	inline vec3 operator + (const vec3& v) const { return vec3(x + v.x, y + v.y, z + v.z); }
	inline vec3 operator - (const vec3& v) const { return vec3(x - v.x, y - v.y, z - v.z); }
	inline vec3 operator * (const vec3& v) const { return vec3(x * v.x, y * v.y, z * v.z); }
	inline vec3 operator / (const vec3& v) const { return vec3(x / v.x, y / v.y, z / v.z); }

	inline vec3 operator + (cfloat s) const { return vec3(x + s, y + s, z + s); }
	inline vec3 operator - (cfloat s) const { return vec3(x - s, y - s, z - s); }
	inline vec3 operator * (float s) const { return vec3(x * s, y * s, z * s); }
	inline vec3 operator / (float s) const { return vec3(x / s, y / s, z / s); }

	inline vec3& operator = (const vec3& v) { x = v.x; y = v.y; z = v.z; return *this; }
	inline vec3& operator += (const vec3& v) { x += v.x; y += v.y; z += v.z; return *this; }
	inline vec3& operator -= (const vec3& v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
	inline vec3& operator *= (const vec3& v) { x *= v.x; y *= v.y; z *= v.z; return *this; }
	inline vec3& operator /= (const vec3& v) { x /= v.x; y /= v.y; z /= v.z; return *this; }

	inline vec3& operator += (cfloat s) { x += s; y += s; z += s; return *this; }
	inline vec3& operator -= (cfloat s) { x -= s; y -= s; z -= s; return *this; }
	inline vec3& operator *= (cfloat s) { x *= s; y *= s; z *= s; return *this; }
	inline vec3& operator /= (cfloat s) { x /= s; y /= s; z /= s; return *this; }

	inline vec3 operator - () const { return vec3(-x, -y, -z); }

	inline float dot2D(const vec3& v) const { return (x * v.x +           z * v.z); }
	inline float dot3D(const vec3& v) const { return (x * v.x + y * v.y + z * v.z); }
	inline float sqLen2D() const { return (dot2D(*this)); }
	inline float sqLen3D() const { return (dot3D(*this)); }
	inline float len2D() const { return sqrtf(sqLen2D()); }
	inline float len3D() const { return sqrtf(sqLen3D()); }

	// cross this vector with <v>
	inline vec3 cross(const vec3& v) const {
		return vec3((y * v.z) - (z * v.y), (z * v.x) - (x * v.z), (x * v.y) - (y * v.x));
	}
	// cross this vector with <v> (in-place)
	inline vec3& icross(const vec3& v) {
		x = (y * v.z) - (z * v.y);
		y = (z * v.x) - (x * v.z);
		z = (x * v.y) - (y * v.x);
		return *this;
	}

	// normalize this vector
	inline vec3 norm() const {
		float dp = sqLen3D();
		float rt = 1.0f;

		if (dp > 0.0f && dp != 1.0f) {
			rt = 1.0f / sqrtf(dp);
		}

		return vec3(x * rt, y * rt, z * rt);
	}

	std::string str() const {
		char buf[64] = {'\0'};
		snprintf(buf, 63, "<%.3f, %.3f, %.3f>", x, y, z);
		// implicit conversion
		return buf;
	}

	float x;
	float y;
	float z;
};

struct vec4: public vec3 {
	vec4(): vec3(0.0f, 0.0f, 0.0f), w(0.0f) {}
	vec4(float _x, float _y, float _z, float _w): vec3(_x, _y, _z), w(_w) {}
	vec4(const vec3& v, float _w): vec3(v), w(_w) {}
	vec4(const vec4& v): vec3(v.x, v.y, v.z), w(v.w) {}

	float w;
};

#endif
