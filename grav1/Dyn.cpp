#include "Dyn.h"

C Dyn::accel(int i, C z) const
{
	C accel;
	for (int j = n - 1; j >= 0; j--)
	{
		if (i == j) continue;
		C zj = tab[j].z, r = zj - z;
		// Newton's gravity has an inherent singularity at r = 0 (zero distance).
		// I have to do something sensible in that situation.
		double w = std::max(abs(r), param.guard0_dist);
		// Newton's law of gravity.
		C aij_newtgrav = param.g * tab[j].m * (1 / w / w / w) * r;
		// Lennard-Jones potential (gradient).
		// Let h = signa(j) / |r|.
		// where j is the "other" particle.
		// Note that h is dimensionless, since both |r| and sigma are lengths.
		// Now, ordinarily, the 8th and 14th powers of h would be calculated,
		// but because of the programming (numerics) issues, specifically
		// the stability of integration over large-enough time slices (dt),
		// I chose different long- and short-range powers.
		double h = tab[j].ljsigma / w,
			h3 = h * h * h, h6 = h3 * h3, h12 = h6 * h6;
		C aij_lj = 4 * param.lj_force_unit * (h * h * h6 - h12) / tab[i].m * r;
		accel += aij_newtgrav + aij_lj;
	}
	return accel;
}

void Dyn::accelall() const
{
	std::unique_ptr<DynEntry[]> copy(new DynEntry[n]{});
	for (int i = n - 1; i >= 0; i--)
	{
		copy[i] = tab[i];
		copy[i].a = accel(i, copy[i].z);
	}
	for (int i = n - 1; i >= 0; i--)
		tab[i] = copy[i];
}

DynEntry Dyn::iter(int i)
{
	// Limit the magnitude of a vector.
	auto limit = [=](double lim, C x) {
		double w = abs(x);
		if (w == 0.0) return x;
		double ww = w / lim,
			mag = lim / w * (2.0 / PI) * atan(ww);
		return mag * x;
		};
	auto speed_limit = [=](C v) { return limit(param.speed_limit, v); };
	auto accel_limit = [=](C a) { return limit(param.accel_limit, a); };

	// Method of leapfrog integration: conserves energy.
	// Unphysical signal filters to prevent instabilities
	// while not doing hard-spheres yet (but I might want to consider doing that).
	double dt = param.dt;
	DynEntry t(tab[i]);
	C v1_(t.v + dt / 2 * t.a), v1(speed_limit(v1_)), z2(t.z + v1 * dt);
	C a2_(accel(i, z2)), a2(accel_limit(a2_));
	C v2_(v1 + dt / 2 * a2), v2(speed_limit(v2_));
	t.z = z2, t.v = v2, t.a = a2;

	// If the particle is too far away from the origin, then replace it.
	// Note: these function pointers are allowed to NOT exist (be null)
	// and thus must be checked for existence before each use.
	if (this->driver.replace && this->driver.test_replace
		&& this->driver.test_replace(t))
	{
		double old_mass = t.m;
		t = this->driver.replace();
		double mass_delta = t.m - old_mass;
		this->mass += mass_delta;
	}

	return t;
}

void Dyn::iterall()
{
	std::unique_ptr<DynEntry[]> copy(new DynEntry[n]{});
	for (int i = n - 1; i >= 0; i--)
		copy[i] = iter(i);
	for (int i = n - 1; i >= 0; i--)
		tab[i] = copy[i];
}

void Dyn::center_vz() const
{
	// Welford's online algorithm.
	// Difficult to parallelize, but at least numerically stable.
	C barycenter, velbias;
	for (int i = n - 1, m = 0; i >= 0; i--)
	{
		double o = ++m;
		barycenter += (tab[i].z * tab[i].m - barycenter) / o;
		if (!std::isfinite(barycenter.real()) || !std::isfinite(barycenter.imag()))
			throw 0;
		velbias += (tab[i].v - velbias) / o;
	}
	barycenter /= mass;
	for (int i = n - 1; i >= 0; i--) tab[i].z -= barycenter, tab[i].v -= velbias;
}
