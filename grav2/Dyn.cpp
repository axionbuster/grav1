#include "Dyn.h"

using std::swap;

void Dyn::precompute()
{
	copy = tab;
	for (int i = n() - 1; i >= 0; i--) m_mass += copy[i].m, copy[i].a = accelerate(i, copy[i].z);
	swap(tab, copy);
}

void Dyn::step()
{
	copy = tab;
	for (int i = n() - 1; i >= 0; i--)
	{
		auto& t = copy[i];
		// Method of leapfrog integration for synchronous velocities.
		// "Kick, drift, kick."
		C const v1 = t.v + par.dt / 2 * t.a, z2 = t.z + par.dt * v1;
		C const a2 = accelerate(i, z2); // N.B.: Can take velocities here, if considered out of step.
		C const v2 = v1 + par.dt / 2 * a2; // N.B.: Enforce synchrony.
		t.z = z2, t.v = v2, t.a = a2;
	}
	swap(tab, copy);
}

void Dyn::bias()
{
	int i, n = this->n();
	C zcm, vcm;
	// Since numerical stability is a problem, use the numerically stable method
	// of Welford's online algorithm to compute the arithmetic mean of the relevant vectors.
#define t tab[i - 1]
#define put(avg, term) avg += ((term) - avg) / (double)i
	for (i = 1; i <= n; i++) put(zcm, t.z * t.m), put(vcm, t.v * t.m);
	zcm /= m_mass, vcm /= m_mass;
	for (i = 1; i <= n; i++) t.z -= zcm, t.v -= vcm;
#undef put
#undef t
}

/// <summary>
/// Compute the acceleration felt by particle at index `i` if it were at
/// the hypothetical location z.
/// </summary>
/// <param name="i">Valid index of the particle</param>
/// <param name="z">The particle's hypothetical location</param>
/// <returns>Acceleration, or force divided by the particle's mass</returns>
C Dyn::accelerate(int i, C const& z) const
{
	if (!drv.pair_force) return 0;
	Entry e(tab[i]); e.z = z;
	C f;
	for (int j = n() - 1; j >= 0; j--) if (i != j) f += drv.pair_force(e, tab[j]);
	return f / e.m;
}
