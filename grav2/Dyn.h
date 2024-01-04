#pragma once
#include "Include.h"
#include <vector>
#include <functional>

class Dyn
{
public:
	struct Param
	{
		double g{};
		double dt{};
	};

	struct Entry
	{
		/// <summary>
		/// Position, velocity, acceleration.
		/// </summary>
		C z, v, a;
		/// <summary>
		/// Mass, radius.
		/// </summary>
		double m{ 1 }, r{ 1 };
	};

	struct Driver
	{
		/// <summary>
		/// Compute the force on the left particle by the right particle.
		/// </summary>
		std::function<C(Entry const&, Entry const&)> pair_force;
	};

	typedef std::vector<Entry> V;

	Dyn(Param const& par, V const& tab)
		: par(par), tab(tab) {}
	Dyn(Param const& par, V&& tab)
		: par(par), tab(std::move(tab)) {}
	Dyn(Param const& par)
		: par(par), tab() {}

	/// <summary>
	/// Exposed set of drivers.
	/// </summary>
	Driver driver;
	/// <summary>
	/// Exposed simulation parameters.
	/// </summary>
	Param par;
	/// <summary>
	/// Exposed dynamical table.
	/// </summary>
	V tab;

	/// <summary>
	/// Precompute all accelerations before the first iteration.
	/// </summary>
	void precompute();

	/// <summary>
	/// Integrate a full time step.
	/// </summary>
	void step();

	/// <summary>
	/// De-bias the positions and velocities
	/// by locating the barycenter at (0, 0) and
	/// setting the average velocity to (0, 0).
	/// </summary>
	void bias();

	/// <summary>
	/// Count the number of particles.
	/// </summary>
	/// <returns></returns>
	int n() const { return tab.size(); }

private:
	/// <summary>
	/// A copy of `tab` made at different times by different
	/// method calls. Only valid within that method call who made it.
	/// </summary>
	V copy;

	/// <summary>
	/// Compute the acceleration felt by particle at index `i` if it were at
	/// the hypothetical location z.
	/// </summary>
	/// <param name="i">Valid index of the particle</param>
	/// <param name="z">The particle's hypothetical location</param>
	/// <returns>Acceleration, or force divided by the particle's mass</returns>
	C accelerate(int i, C const& z) const
	{
		if (!driver.pair_force) throw std::exception("pair_force() must exist");
		Entry e(tab[i]); e.z = z;
		C f;
		int n = this->n();
		for (int j = n - 1; j >= 0; j--)
		{
			if (i == j) continue;
			f += driver.pair_force(e, tab[j]);
		}
		return f / e.m;
	}
};

