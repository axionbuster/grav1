#include "Dyn.h"
#include "Beasons.h"

using std::swap;

// Names (usually):
// e: Entry (particle).
//    A particle is called an "entry" because it's a row ("entry") in the dynamical table (`tab`)
//    which is like a spreadsheet of all relevant properties describing each particle.
// tab: Dynamical table.
//    Explained above.
// ...cm: ... of the center of mass or momentum.
// f: Force.

void Dyn::precompute()
{
	copy = tab;
	for (int i = n() - 1; i >= 0; i--)
	{
		auto& e = copy[i];
		m_mass += e.m;
		m_area += e.r * e.r * PI64;
		e.a = accelerate(i);
	}
	swap(tab, copy);
}

void Dyn::step()
{
	// Choice of integrator.
	using namespace beasons;

	// Count cases that are too inaccurate.
	int go_finer{};

	// Count cases that suggest integration step size may be safely increased.
	int go_coarser{};

	copy = tab;
	for (int i = n() - 1; i >= 0; i--)
	{
		auto& e = copy[i];

		// ::: Beason's method of integration with step size adjustment. :::

		// For this entry only, try with this step size.
		// Forgotten at end of each iteration of the for-loop here.
		double dt = par.dt;
		// Number of times to retry in the worst case.
		// If 0, give up.
		int motivation = 4;
		// Latest results from the integrator.
		beasons::BeasonsResults aa;

		do
		{
			// 1. Do the math.

			auto accel = [&](C z, C v)
				{
					e.z = z, e.v = v; // Destructive modification of the copied entry.
					return accelerate(i); // e and i refer to the same entry.
				};
			aa = beasons::beason_bogacki_shampine(par.dt, accel, e.z, e.v, e.a);

			// 2. Quality control (adjust step size).

			// Map signed integers to {-1, 0, 1}, according to the sign.
			// Note 1: (!) is boolean negation. (!!) is doing that twice,
			// which maps all integers to {0, 1}: the output is 0 iff
			// the number is 0.
			// Note 2: As it happens in math, (+1) - 2 = (-1).
			auto sign = [](int a) { return -2 * !!(a < 0) + !!a; };

			if (drv.judge_z)
			{
				// after `sign` call:
				// -1: inhibition, try finer time step.
				// 0: neutral, do nothing in particular.
				// +1: ambition, try coarser time step.
				switch (sign(drv.judge_z(aa.y0_strong, aa.y0_weak)))
				{
				case 0: // pass-through
				case -1:
					dt = std::max(par.low_dt, dt / 2);
					go_finer++;
					// Retry with less "motivation."
					// Too little motivation causes the loop to just give up.
					motivation--;
					continue;
				case 1:
					go_coarser++;
					// No need to adjust `dt` since we're moving onto
					// another particle.
					break;
				}
			}
			if (drv.judge_v)
			{
				switch (sign(drv.judge_v(aa.y1_strong, aa.y1_weak)))
				{
				case 0: // pass-through
				case -1:
					dt = std::max(par.low_dt, dt / 2);
					go_finer++, motivation--;
					continue;
				case 1:
					go_coarser++;
					break;
				}
			}
			break;
		} while (motivation);

		e.z = aa.y0_strong;
		e.v = aa.y1_strong;
		e.a = aa.y2;
	}

	// Apply *global* time step adjustment.
	if (go_finer) par.dt = std::max(par.low_dt, par.dt / 2);
	else if (go_coarser) par.dt = std::min(par.high_dt, par.dt * 2);

	swap(tab, copy);
}

void Dyn::bias()
{
	// Barycenter and momentum.
	C zcm, vcm;
	// Since numerical stability is a problem, use the numerically stable method
	// of Welford's online algorithm to compute the arithmetic mean of the relevant vectors.
#define e tab[i - 1]
#define put(avg, term) avg += ((term) - avg) / (double)i
	for (int i = 1; i <= n(); i++) put(zcm, e.z * e.m), put(vcm, e.v * e.m);
	zcm /= m_mass, vcm /= m_mass;
	for (int i = 1; i <= n(); i++) e.z -= zcm, e.v -= vcm;
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
