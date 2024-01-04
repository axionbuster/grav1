#include "Dyn.h"

void Dyn::precompute()
{
}

void Dyn::step()
{
}

void Dyn::bias()
{
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
	int n = this->n(); C f;
	for (int j = n - 1; j >= 0; j--)
		if (i != j) f += drv.pair_force(e, tab[j]);
	return f / e.m;
}
