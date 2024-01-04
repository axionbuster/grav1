#include "Dyn.h"

using std::swap;

void Dyn::precompute()
{
	copy = tab;
	for (int i = n() - 1; i >= 0; i--) m_mass += copy[i].m, copy[i].a = accelerate(i);
	swap(tab, copy);
}

void Dyn::step()
{
	copy = tab;
	for (int i = n() - 1; i >= 0; i--)
	{
		auto& e = copy[i];
		// Method of leapfrog integration for synchronous velocities.
		// "Kick, drift, kick."
		C const v1 = e.v + par.dt / 2 * e.a;
		e.z += par.dt * v1;
		e.a = accelerate(i); // N.B.: Can take velocities here, if considered out of step.
		e.v = v1 + par.dt / 2 * e.a; // N.B.: Enforce synchrony.
	}
	swap(tab, copy);
}

void Dyn::bias()
{
	int i, n = this->n();
	C zcm, vcm;
	// Since numerical stability is a problem, use the numerically stable method
	// of Welford's online algorithm to compute the arithmetic mean of the relevant vectors.
#define e tab[i - 1]
#define put(avg, term) avg += ((term) - avg) / (double)i
	for (i = 1; i <= n; i++) put(zcm, e.z * e.m), put(vcm, e.v * e.m);
	zcm /= m_mass, vcm /= m_mass;
	for (i = 1; i <= n; i++) e.z -= zcm, e.v -= vcm;
#undef put
#undef e
}

/// <summary>
/// Compute the acceleration felt by particle at index `i` if it were at
/// the hypothetical location z.
/// </summary>
/// <param name="i">Valid index of the particle</param>
/// <param name="z">The particle's hypothetical location</param>
/// <returns>Acceleration, or force divided by the particle's mass</returns>
C Dyn::accelerate(int i) const
{
	if (!drv.pair_force) return 0;
	Entry e(tab[i]);
	C f;
	for (int j = n() - 1; j >= 0; j--) if (i != j) f += drv.pair_force(e, tab[j]);
	return f / e.m;
}
