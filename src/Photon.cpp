#include <cmath>
#include <cassert>

#include "Photon.hpp"
#include "vec3.hpp"
#include "matrix44.hpp"
#include "RNG.hpp"
#include "Dust.hpp"
#include "Slab.hpp"
#include "Grid.hpp"
#include "Image.hpp"
#include "Observer.hpp"

Photon::Photon(const RNG& rng, const vec3& p, const Grid& g):
	pos(p.x,  p.y,  p.z),
	dir(0.0f, 0.0f, 1.0f),
	pol(1.0f, 0.0f, 0.0f, 0.0f),
	numScatEvents(0), absorbed(false)
{
	cfloat xi0 = rng.FUniformRandom();
	cfloat xi1 = rng.FUniformRandom();

	UpdateAngles(acosf(sqrtf(xi0)), M_PI2 * xi1);
	UpdateDir();
	UpdateCellIndices(g);
}

void Photon::Scatter(const RNG& rng, const Dust& dust, ScatterMode sm) {
	numScatEvents++;

	switch (sm) {
		case SM_ISOTROPIC: { IsotropicScatter(rng, dust); } break;
		// case SM_ANISOTROPIC: { AnisotropicScatter(rng, dust); } break;
		case SM_ANISOTROPIC: { AnisotropicScatterBW(rng, dust); } break;
	}

	UpdateDir();

	#ifdef DEBUG
	if (isnan( the)) { assert(false); }
	if (isnan( phi)) { assert(false); }
	if (isnan(cthe)) { assert(false); }
	if (isnan(sthe)) { assert(false); }
	if (isnan(cphi)) { assert(false); }
	if (isnan(sphi)) { assert(false); }

	if (isnan(pol.x)) { assert(false); }
	if (isnan(pol.y)) { assert(false); }
	if (isnan(pol.z)) { assert(false); }
	if (isnan(pol.w)) { assert(false); }
	#endif
}

void Photon::IsotropicScatter(const RNG& rng, const Dust&) {
	cfloat xi0 = rng.FUniformRandom();
	cfloat xi1 = rng.FUniformRandom();

	UpdateAngles(acosf((2.0f * xi0) - 1.0f), M_PI2 * xi1);
}

void Photon::AnisotropicScatter(const RNG& rng, const Dust& dust) {
	cfloat g  = dust.g;
	cfloat gg = dust.gg;
	cfloat pl = dust.pl;
	cfloat pc = dust.pc;
	cfloat sc = dust.sc;

	// Henyey-Greenstein phase function
	// equation 11 from "The DIRTY Model"
	// equation 33 from "Intro. to MC RT"
	//
	// NOTE: if g is 0, this is DIVZERO
	// so we cannot simulate fully iso-
	// tropic scattering
	cfloat n = 1.0f / (2.0f * g);
	cfloat d =
		(1.0f - gg) /
		(1.0f - g + 2.0f * g * rng.FUniformRandom());
	cfloat cTHE = std::max(-1.0f, std::min(1.0f, n * ((1.0f + gg) - (d * d))));

	// cos(T) = HenyeyGreenstein() ==>
	// T = arccos(HenyeyGreenstein())
	cfloat  THE  = acosf(cTHE);
	cfloat sTHE  = sinf(THE);
	cfloat  THEf = THE * (1.0f + 3.13f * sc * expf((-7.0f * THE) / M_PI));
	cfloat cTHEf = cosf(THEf);

	#ifdef DEBUG
	assert(!isnan(cTHE));
	assert(!isnan(sTHE));
	assert(!isnan( THE));
	#endif



	// set the new phi and its *sines
	cfloat nphi = M_PI2 * rng.FUniformRandom();
	cfloat dphi = phi - nphi;

	phi = nphi;
	cphi = cosf(phi);
	sphi = sinf(phi);

	cfloat tx = (dir.x * cphi * cthe - dir.y * sphi) / (sthe != 0.0f? sthe: EPSILON);
	cfloat ty = (dir.y * cphi * cthe + dir.x * sphi) / (sthe != 0.0f? sthe: EPSILON);
	cfloat tz = sthe;
	vec3 t(tx, ty, tz);
	// set the new Cartesian direction vector
	// NOTE: will be overwritten by UpdateDir()
	dir = (dir * cTHE + t * sTHE).norm();


	cuint sign = (M_PI < dphi && dphi < M_PI2)? 1: -1;
	cfloat i1div = sign * (sqrtf(1.0f - cTHE * cTHE) * sqrtf(1.0f - dir.z * dir.z));
	cfloat i2div = sign * (sqrtf(1.0f - cTHE * cTHE) * sqrtf(1.0f -  cthe *  cthe));
	cfloat cosi1 = (-cthe  + dir.z * cTHE) / (i1div != 0.0f? i1div: EPSILON);
	cfloat cosi2 = (-dir.z + cthe  * cTHE) / (i2div != 0.0f? i2div: EPSILON);
	cfloat i1 = acosf(std::max(-1.0f, std::min(cosi1, 1.0f)));
	cfloat i2 = acosf(std::max(-1.0f, std::min(cosi2, 1.0f)));
	cfloat ci1 = cosf(2.0f * -i1), ci2 = cosf(2.0f * (M_PI - i2));
	cfloat si1 = sinf(2.0f * -i1), si2 = sinf(2.0f * (M_PI - i2));

	#ifdef DEBUG
	assert(!isnan(i1));
	assert(!isnan(i2));
	assert(!isnan(ci1) && !isnan(ci2));
	assert(!isnan(si1) && !isnan(si2));
	#endif

	// set the new theta and its *sines
	the = acosf(dir.z);
	cthe = dir.z;
	sthe = sinf(the);



	// with T (and cos(T)) now known, set the R-matrix parameters
	// (these describe the scattering probability in the frame of
	// the dust particle with respect to the incident direction)
	cfloat P1 = ((1.0f - gg) / powf((1.0f + gg - 2.0f * g * cTHE), 1.5f));
	cfloat P2 = (-pl * P1 * ((1.0f - cTHE * cTHE) / (1.0f + cTHE * cTHE)));
	cfloat P3 = (P1 * ((2.0f * cTHE) / (1.0f + cTHE * cTHE)));
	cfloat P4 = (-pc * P1 * ((1.0f - cTHEf * cTHEf) / (1.0f + cTHEf * cTHEf)));

	#ifdef DEBUG
	assert(!isnan(P1));
	assert(!isnan(P2));
	assert(!isnan(P3));
	assert(!isnan(P4));
	#endif

	// compose R(T) for spherical particles
	matrix44 R;
	R.setCol(0, vec4(  P1,   P2, 0.0f, 0.0f)); // 1st column
	R.setCol(1, vec4(  P2,   P1, 0.0f, 0.0f)); // 2nd column
	R.setCol(2, vec4(0.0f, 0.0f,   P3,   P4)); // 3rd column
	R.setCol(3, vec4(0.0f, 0.0f,  -P4,   P3)); // 4th column
	R.mul(0.75f);





	// NOTE: should i1 REALLY be a random number?
	// "The DUSTY Model" mentions that "angles i1
	// and i2 are calculated from the knowledge of
	// the direction the photon is traveling and
	// theta(obs) using spherical geometry", but
	// "Intro. to MC RT" says "i1 is sampled from
	// a uniform distribution 2PI * Random(0, 1)"?


	// compose the matrix L(-i1)
	matrix44 L1;
	L1.setCol(0, vec4(1.0f, 0.0f, 0.0f, 0.0f)); // 1st column (a, e, i, m)
	L1.setCol(1, vec4(0.0f,  ci1, -si1, 0.0f)); // 2nd column (b, f, j, n)
	L1.setCol(2, vec4(0.0f,  si1,  ci1, 0.0f)); // 3rd column (c, g, k, o)
	L1.setCol(3, vec4(0.0f, 0.0f, 0.0f, 1.0f)); // 4th column (d, h, l, p)
	// L1.mul(0.75f);

	// compose the matrix L(PI - i2)
	matrix44 L2;
	L2.setCol(0, vec4(1.0f, 0.0f, 0.0f, 0.0f)); // 1st column (a, e, i, m)
	L2.setCol(1, vec4(0.0f,  ci2, -si2, 0.0f)); // 2nd column (b, f, j, n)
	L2.setCol(2, vec4(0.0f,  si2,  ci2, 0.0f)); // 3rd column (c, g, k, o)
	L2.setCol(3, vec4(0.0f, 0.0f, 0.0f, 1.0f)); // 4th column (d, h, l, p)
	// L2.mul(0.75f);

	// compute S' = L(i2) * R(T) * L(-i1) * S
	pol = L2.mul(R.mul(L1.mul(pol)));

	// "the Stokes polarization vector is
	// renormalized after each scattering
	// so that I = 1" ("The DUSTY Model")
	// pol.y /= pol.x;
	// pol.z /= pol.x;
	// pol.w /= pol.x;
	// pol.x  = 1.0f;
	//
	// "the scattered Stokes parameters are
	// normalized to P1" ("Intro. to MC RT")
	pol.x /= P1;
	pol.y /= P1;
	pol.z /= P1;
	pol.w /= P1;
}



// taken from the example Fortran programs by Bjorkman et all
// FIXME: effect of changing g seems very wrong, no rings etc.
void Photon::AnisotropicScatterBW(const RNG& rng, const Dust& dust) {
	cfloat g = dust.g, gg = dust.gg;

	float tmp, p1, p2, p3, p4;
	float a11, a12, a13, a21, a22, a23, a24, a31, a32, a33, a34, a42, a43, a44;
	float ri3, cos2, sin2, phip;
	float bott, cosi1, cosi2, cosi3, sini1, sini2, sini3, cos2i1, cos2i2, cos2i3, sin2i1, sin2i2, sin2i3, costp, sintp, sin2cos1, cos2sin1;

	// normalize the Stokes vector
	cfloat w = pol.x;
	pol.x /= w;
	pol.y /= w;
	pol.z /= w;
	pol.w /= w;

	// copy the initial variables
	costp = cthe;
	sintp = sthe;
	phip = phi;

	// cos(THETA) according to (33)
	cfloat n = 1.0f / (2.0f * g);
	cfloat d =
		(1.0f - gg) /
		(1.0f - g + 2.0f * g * rng.FUniformRandom());
	// prevent NaN's in acosf()
	cfloat cTHE = std::max(-1.0f, std::min(1.0f, n * ((1.0f + gg) - (d * d))));

	cfloat cTHEsq = cTHE * cTHE;
	cfloat sTHE = sqrtf(1.0f - cTHEsq);

	// If the scattering angle is 0, return,
	// if it is PI just invert the U element
	if (cTHE == 1.0f || cTHE == -1.0f) {
		if (cTHE == -1.0f) {
			// NOTE: is this enough? need to
			// adjust cthe and sthe at least?
			pol.z = -pol.z;
		}
		return;
	}

	// calculate P1, ..., P4 of the dust matrix
	dust.GetMatrixParameters(p1, p2, p3, p4, cTHE, cTHEsq);

	cfloat ri1 = M_PI2 * rng.FUniformRandom();


	// TODO: refactor
	if (ri1 > M_PI) {
		// reduce ri1 to [0, pi]
		ri3 = M_PI2 - ri1;

		cosi3 = cosf(ri3);
		sini3 = sinf(ri3);
		sin2i3 = sini3 * 2.0f * cosi3;
		cos2i3 = cosi3 * 2.0f * cosi3 - 1.0f;

		// (cosine of THETA) in global coordinates;
		// THETA is the photon-local zenith SA
		cthe = costp * cTHE + sintp * sTHE * cosi3;

		if (fabs(cthe) < 1.0f) {
			tmp = sqrtf(1.0f - cthe * cthe);
			sthe = fabs(tmp);
			sini2 = sini3 * sintp / sthe; // (sini1 * sintp / sthe) in else-branch
			bott = sthe * sTHE;
			cosi2 = costp / bott - cthe * cTHE / bott;
		} else {
			sthe = 0.0f;
			sini2 = 0.0f;
			if (cthe >=  1.0f) { cosi2 = -1.0f; }
			if (cthe <= -1.0f) { cosi2 =  1.0f; }
		}

		// (cosine of PHI) in global coordinates;
		// PHI is the photon-local azimuthal SA
		float cPHI = -cosi2 * cosi3 + sini2 * sini3 * cTHE;
		if (cPHI >  1.0f) { cPHI =  1.0f; }
		if (cPHI < -1.0f) { cPHI = -1.0f; }

		UpdateAngles(acosf(cthe), phip + acosf(cPHI));

		sin2i2 = sini2 * 2.0f * cosi2;
		cos2i2 = cosi2 * 2.0f * cosi2 - 1.0f;
		sin2 = sin2i2 * sin2i3;
		cos2 = cos2i2 * cos2i3;
		sin2cos1 = sin2i2 * cos2i3;
		cos2sin1 = cos2i2 * sin2i3;

		// elements of the composite matrix L(pi - i2)R(T)L(-i1)
		a11 = p1;           a12 = p2 * cos2i3;                    a13 = p2 * sin2i3;
		a21 = p2 * cos2i2;  a22 = p1 * cos2 - p3 * sin2;          a23 = p1 * cos2sin1 + p3 * sin2cos1; a24 = -p4 * sin2i2;
		a31 = -p2 * sin2i2; a32 = -p1 * sin2cos1 - p3 * cos2sin1; a33 = -p1 * sin2 + p3 * cos2;        a34 = -p4 * cos2i2;
		                    a42 = -p4 * sin2i3;                   a43 = p4 * cos2i3;                   a44 = p3;
	} else {
		cosi1 = cosf(ri1);
		sini1 = sinf(ri1);
		sin2i1 = sini1 * 2.0f * cosi1;
		cos2i1 = cosi1 * 2.0f * cosi1 - 1.0f;

		cthe = costp * cTHE + sintp * sTHE * cosi1;

		if (fabs(cthe) < 1.0f) {
			tmp = sqrtf(1.0f - cthe * cthe);
			sthe = fabs(tmp);
			sini2 = sini1 * sintp / sthe; // (sini3 * sintp / sthe) in if-branch
			bott = sthe * sTHE;
			cosi2 = costp / bott - cthe * cTHE / bott;
		} else {
			sthe = 0.0f;
			sini2 = 0.0f;
			if (cthe >=  1.0f) { cosi2 = -1.0f; }
			if (cthe <= -1.0f) { cosi2 =  1.0f; }
		}

		// (cosine of PHI) in global coordinates;
		// PHI is the photon-local azimuthal SA
		float cPHI = -cosi1 * cosi2 + sini1 * sini2 * cTHE;
		if (cPHI >  1.0f) { cPHI =  1.0f; }
		if (cPHI < -1.0f) { cPHI = -1.0f; }

		UpdateAngles(acosf(cthe), phip + acosf(cPHI));

		sin2i2 = sini2 * 2.0f * cosi2;
		cos2i2 = cosi2 * 2.0f * cosi2 - 1.0f;
		sin2 = sin2i2 * sin2i1;
		cos2 = cos2i2 * cos2i1;
		sin2cos1 = sin2i2 * cos2i1;
		cos2sin1 = cos2i2 * sin2i1;

		a11 = p1;          a12 = p2 * cos2i1;                   a13 = -p2 * sin2i1;
		a21 = p2 * cos2i2; a22 = p1 * cos2 - p3 * sin2;         a23 = -p1 * cos2sin1 - p3 * sin2cos1; a24 = p4 * sin2i2;
		a31 = p2 * sin2i2; a32 = p1 * sin2cos1 + p3 * cos2sin1; a33 = -p1 * sin2 + p3 * cos2;         a34 = -p4 * cos2i2;
		                   a42 = p4 * sin2i1;                   a43 = p4 * cos2i1;                    a44 = p3;
	}

	// calculate (L2 * R * L1) * [I, Q, U, V]
	pol.x = (a11 * pol.x + a12 * pol.y + a13 * pol.z) / p1;
	pol.y = (a21 * pol.x + a22 * pol.y + a23 * pol.z + a24 * pol.w) / p1;
	pol.z = (a31 * pol.x + a32 * pol.y + a33 * pol.z + a34 * pol.w) / p1;
	pol.w = (a42 * pol.y + a43 * pol.z + a44 * pol.w) / p1;

	// set the final Stokes vector
	pol.x *= w;
	pol.y *= w;
	pol.z *= w;
	pol.w *= w;
}



// theta and phi are assumed to be in radians
void Photon::UpdateAngles(cfloat t, cfloat p) {
	the = t;
	phi = p;

	if (the <  0.0f) { the += M_PI2; }
	if (the > M_PI2) { the -= M_PI2; }
	if (phi <  0.0f) { phi += M_PI2; }
	if (phi > M_PI2) { phi -= M_PI2; }

	cthe = cosf(the); sthe = sinf(the);
	cphi = cosf(phi); sphi = sinf(phi);

	// #ifdef DEBUG
	// assert(!isnan(the) && !isnan(phi));
	// #endif
}

// converts angles in spherical coordinates to
// a normalized 3D Cartesian direction vector
// (a point on the surface of the unit sphere)
//
// angles are defined with respect to the "OGL"
// coordinate system: theta is the zenith angle
// relative to y+ (and rotates clockwise around
// z if phi is 0), phi is azimuth angle relative
// to x+ (and rotates clockwise around y if theta
// is 0)
//
//           y+
//         |
//         |
//         |
//         |__________ x+
//        /
//      /  
//    /
//   z+
//
void Photon::UpdateDir() {
	// Woods et all (implicitly defines theta
	// wrt. z-axis and phi wrt. x- or y-axis)
	dir.x = 1.0f * sthe * cphi;
	dir.y = 1.0f * sthe * sphi;
	dir.z = 1.0f * cthe;

	// this breaks UpdateIntensityMoments()
	// dir.x = 1.0f * sthe * cphi;
	// dir.y = 1.0f * cthe;
	// dir.z = 1.0f * sthe * sphi;
}

void Photon::UpdatePos(cfloat d) {
	pos += (dir * d);
}

void Photon::UpdateCellIndices(const Grid& g) {
	xcell = g.GetXCellIdx(pos.x);
	ycell = g.GetYCellIdx(pos.y);
	zcell = g.GetZCellIdx(pos.z);

	#ifdef DEBUG
	assert(g.IndicesInBounds(xcell, ycell, zcell));
	#endif
}

std::string Photon::str() {
	char buf1[512] = {'\0'}; snprintf(buf1, 511, "the: %.3f, cthe: %.3f, sthe: %.3f", the, cthe, sthe);
	char buf2[512] = {'\0'}; snprintf(buf2, 511, "phi: %.3f, cphi: %.3f, sphi: %.3f", phi, cphi, sphi);
	return (dir.str() + "\t" + std::string(buf1) + "\t" + std::string(buf2));
}



bool Photon::PerspProject(const Observer& o, Image* img) const {
	// plug [P(t) = P(0) + D * t] into plane eq.
	cfloat A = (-o.zdir).dot3D(pos);
	cfloat B = (-o.zdir).dot3D(dir);
	cfloat t = ((o.dst - A) / B);

	if (t < 0.0f) {
		return false;
	}

	// photon ray-plane intersection pos
	const vec3 p = pos + (dir * t);
	const uint ix = uint( p.x + img->hw);
	const uint iy = uint(-p.y + img->hh);

	if (ix < 1 || ix >= img->w - 1) { return false; }
	if (iy < 1 || iy >= img->h - 1) { return false; }

	// "over-brighten" pixels in a 3x3 neighbourhood to counter
	// under-sampling (disabled since not of any scientific use)
	// img->d[ix - 1][iy - 1] = std::min(PXL_MAX_VAL, img->d[ix - 1][iy - 1] + 1);   img->m = std::max(img->m, img->d[ix - 1][iy - 1]);
	// img->d[ix    ][iy - 1] = std::min(PXL_MAX_VAL, img->d[ix    ][iy - 1] + 2);   img->m = std::max(img->m, img->d[ix    ][iy - 1]);
	// img->d[ix + 1][iy - 1] = std::min(PXL_MAX_VAL, img->d[ix + 1][iy - 1] + 1);   img->m = std::max(img->m, img->d[ix + 1][iy - 1]);
	// img->d[ix - 1][iy    ] = std::min(PXL_MAX_VAL, img->d[ix - 1][iy    ] + 2);   img->m = std::max(img->m, img->d[ix - 1][iy    ]);
	img->d[ix    ][iy    ] = std::min(PXL_MAX_VAL, img->d[ix    ][iy    ] + 4);   img->m = std::max(img->m, img->d[ix    ][iy    ]);
	// img->d[ix + 1][iy    ] = std::min(PXL_MAX_VAL, img->d[ix + 1][iy    ] + 2);   img->m = std::max(img->m, img->d[ix + 1][iy    ]);
	// img->d[ix - 1][iy + 1] = std::min(PXL_MAX_VAL, img->d[ix - 1][iy + 1] + 1);   img->m = std::max(img->m, img->d[ix - 1][iy + 1]);
	// img->d[ix    ][iy + 1] = std::min(PXL_MAX_VAL, img->d[ix    ][iy + 1] + 2);   img->m = std::max(img->m, img->d[ix    ][iy + 1]);
	// img->d[ix + 1][iy + 1] = std::min(PXL_MAX_VAL, img->d[ix + 1][iy + 1] + 1);   img->m = std::max(img->m, img->d[ix + 1][iy + 1]);

	return true;

	// cfloat rx = (2.0f * o.dst) / img->w;
	// cfloat ry = (2.0f * o.dst) / img->h;
	// cuint ix = uint((( pos.x + o.dst) / rx) + img->hw);
	// cuint iy = uint(((-pos.y + o.dst) / ry) + img->hh);
}

bool Photon::OrthoProject(const Observer& o, Image* img) const {
	if (dir.dot3D(o.zdir) < -0.95f) {
		// only project photons which exited more or
		// less parallel to the viewing direction, no
		// perspective division (note that this only
		// produces useful images for large grids in
		// which photons get the chance to spread out
		// prior to projection)
		cuint ix =  pos.x + img->hw;
		cuint iy = -pos.y + img->hh;

		if (ix < 1 || ix >= img->w - 1) { return false; }
		if (iy < 1 || iy >= img->h - 1) { return false; }

		// img->d[ix - 1][iy - 1] = std::min(PXL_MAX_VAL, img->d[ix - 1][iy - 1] + 1);   img->m = std::max(img->m, img->d[ix - 1][iy - 1]);
		// img->d[ix    ][iy - 1] = std::min(PXL_MAX_VAL, img->d[ix    ][iy - 1] + 2);   img->m = std::max(img->m, img->d[ix    ][iy - 1]);
		// img->d[ix + 1][iy - 1] = std::min(PXL_MAX_VAL, img->d[ix + 1][iy - 1] + 1);   img->m = std::max(img->m, img->d[ix + 1][iy - 1]);
		// img->d[ix - 1][iy    ] = std::min(PXL_MAX_VAL, img->d[ix - 1][iy    ] + 2);   img->m = std::max(img->m, img->d[ix - 1][iy    ]);
		img->d[ix    ][iy    ] = std::min(PXL_MAX_VAL, img->d[ix    ][iy    ] + 4);   img->m = std::max(img->m, img->d[ix    ][iy    ]);
		// img->d[ix + 1][iy    ] = std::min(PXL_MAX_VAL, img->d[ix + 1][iy    ] + 2);   img->m = std::max(img->m, img->d[ix + 1][iy    ]);
		// img->d[ix - 1][iy + 1] = std::min(PXL_MAX_VAL, img->d[ix - 1][iy + 1] + 1);   img->m = std::max(img->m, img->d[ix - 1][iy + 1]);
		// img->d[ix    ][iy + 1] = std::min(PXL_MAX_VAL, img->d[ix    ][iy + 1] + 2);   img->m = std::max(img->m, img->d[ix    ][iy + 1]);
		// img->d[ix + 1][iy + 1] = std::min(PXL_MAX_VAL, img->d[ix + 1][iy + 1] + 1);   img->m = std::max(img->m, img->d[ix + 1][iy + 1]);

		return true;
	}

	return false;
}
