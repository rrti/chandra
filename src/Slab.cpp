#include "RNG.hpp"
#include "Slab.hpp"
#include "Photon.hpp"
#include "ThreadData.hpp"

Slab::Slab(cfloat _zmin, cfloat _zmax, cfloat _den, cfloat _opa, cuint _numLvls):
	zmin(_zmin), zmax(_zmax),
	zRange(zmax - zmin),
	density(_den),
	opacity(_opa),
	tauMax(density * opacity),
	numLevels(_numLvls) {
}

// note: this implements only the summation parts of J, H, and K
void Slab::UpdateIntensityMoments(ThreadData* td, const vec3& opos, const Photon& p) {
	cfloat z1 =  opos.z;
	cfloat z2 = p.pos.z;
	uint l1 = 0, l2 = 0;

	if ((z1 > zmin && z2 > zmin) && (uint(z1 / zRange * numLevels) == uint(z2 / zRange * numLevels))) {
		// the photon did not cross over to another
		// slab slice on its last position update
		return;
	}

	// note: this might fail if zmin < 0.0f
	// since (zn / zRange) would be negative
	if (p.cthe > 0.0f) {
		// update the positive moments
		if (z1 <= zmin) {
			l1 = 0;
		} else {
			l1 = uint((z1 / zRange) * (numLevels));
		}

		if (z2 >= zRange) {
			l2 = numLevels - 1;
		} else {
			l2 = uint((z2 / zRange) * (numLevels));
		}

		// #ifdef DEBUG
		// assert(l2 < numLevels);
		// #endif

		for (uint i = l1; i <= l2; i++) {
			// (cthe / abs(cthe)) == 1 so in this case
			// hPos[i] += 1; similarly kPos[i] += cthe
			// since (cthe^2 / abs(cthe)) == cthe
			td->jPos[i] += (1.0f / p.cthe);
			td->hPos[i] += 1.0f;
			td->kPos[i] += p.cthe;
		}
	} else {
		// update the negative moments
		if (z1 >= zRange) {
			l1 = numLevels - 1;
		} else {
			l1 = uint((z1 / zRange) * (numLevels));
		}

		if (z2 <= zmin) {
			l2 = 0;
		} else {
			l2 = uint((z2 / zRange) * (numLevels));
		}

		// #ifdef DEBUG
		// assert(l1 < numLevels);
		// #endif

		for (uint i = l2; i <= l1; i++) {
			td->jNeg[i] += (1.0f / -p.cthe);
			td->hNeg[i] -= 1.0f;
			td->kNeg[i] += -p.cthe;
		}
	}
}

float Slab::GetRandomL(const RNG& rng) const {
	// sample a random optical depth
	// and non-interaction length d
	cfloat t = -logf(rng.FUniformRandom() + 0.0001f);
	cfloat d = (t * zRange) / tauMax;
	return d;
}
