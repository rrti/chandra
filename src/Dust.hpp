#ifndef DUST_HDR
#define DUST_HDR

#include "TypeDefs.hpp"

struct Dust {
	Dust(cfloat _a, cfloat _g, cfloat _pl, cfloat _pc, cfloat _sc):
		a(_a), g(_g), gg(g * g), pl(_pl), pc(_pc), sc(_sc) {
	}

	void GetMatrixParameters(float& p1, float& p2, float& p3, float& p4, cfloat cTHE, cfloat cTHEsq) const {
		cfloat tmp = (gg + 1.0f - g * 2.0f * cTHE);
		p1 = (1.0f - gg) / powf(tmp, 1.5f);
		p2 = -(pl) * p1 * (1.0f - cTHEsq) / (cTHEsq + 1.0f);
		p3 = p1 * 2.0f * cTHE / (cTHEsq + 1.0f);

		// get THETA in degrees
		cfloat THE = acosf(cTHE) * RAD2DEG;
		cfloat THEf_deg = THE * 3.13f * expf(THE * -7.0f / 180.0f);
		cfloat THEf_rad = (THE + sc * THEf_deg) * DEG2RAD;

		cfloat cTHEf = cos(THEf_rad);
		cfloat cTHEf_sq = cTHEf * cTHEf;

		p4 = -(pc) * p1 * (1.0f - cTHEf_sq) / (cTHEf_sq + 1.0f);
	}

	// if we are considering a purely scattering
	// atmosphere, albedo should be 1 (making all
	// interactions scatter-events)
	//
	// note: all these parameters are assumed to
	// be equal across all cells in the 3D grid
	// (only density and opacity can vary), ie.
	// we have "the same type of dust everywhere"
	float a;		// albedo
	float g;		// scattering asymmetry parameter g = <cos(theta)>
	float gg;		// squared value of g
	float pl;		// peak linear polarization
	float pc;		// peak circular polarization
	float sc;		// skew factor
};

#endif
