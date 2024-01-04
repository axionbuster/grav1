#pragma once
#include "Include.h"
#include <vector>
#include <functional>

/// <summary>
/// Simulation.
/// </summary>
class Dyn
{
public:
	/// <summary>
	/// Simulation parameters in world units.
	/// </summary>
	struct Param
	{
		/// <summary>
		/// Universal gravitational constant (units: LLL/T/T/M).
		/// </summary>
		double g{};
		/// <summary>
		/// Time step (T per step).
		/// </summary>
		double dt{};
	};

	/// <summary>
	/// Kinematic and dynamic properties.
	/// </summary>
	struct Entry
	{
		/// <summary>
		/// Position (L), velocity (L/T), acceleration (L/T/T).
		/// </summary>
		C z, v, a;
		/// <summary>
		/// Mass (M), radius (L).
		/// </summary>
		double m{ 1 }, r{ 1 };
	};

	/// <summary>
	/// Externally specified behavior.
	/// </summary>
	struct Driver
	{
		/// <summary>
		/// Compute the force on the left (first) particle by the right (second) particle.
		/// </summary>
		std::function<C(Param const &par, Entry const&, Entry const&)> pair_force;
	};

	typedef std::vector<Entry> V;

	// For these:
	// Anything except `tab` doesn't make sense to be "moved" without being copied.

	Dyn() = default;
	Dyn(Param const& par) : par(par) {}
	Dyn(Dyn const& dyn) = default;
	Dyn(Dyn&& dyn) noexcept : par(dyn.par), tab(std::move(dyn.tab)), drv(dyn.drv), copy(), m_mass(dyn.m_mass) {}
	Dyn& operator=(Dyn const& dyn) noexcept
	{
		if (&dyn == this) return *this;
		par = dyn.par, tab = dyn.tab, drv = dyn.drv, m_mass = dyn.m_mass;
		copy = V();
		return *this;
	}
	Dyn& operator=(Dyn&& dyn) noexcept
	{
		if (&dyn == this) return *this;
		par = dyn.par, drv = dyn.drv, m_mass = dyn.m_mass;
		tab = std::move(dyn.tab);
		copy = V();
		return *this;
	}
	Entry& operator[](int i) { return tab[i]; }
	Entry const& operator[](int i) const { return tab[i]; }

	/// <summary>
	/// Exposed simulation parameters.
	/// </summary>
	Param par;
	/// <summary>
	/// Exposed dynamical table. (Stores all kinematical and dynamical variables
	/// of all particles). Do not exceed the size of `int` (signed).
	/// </summary>
	V tab;
	/// <summary>
	/// Exposed set of drivers.
	/// </summary>
	Driver drv;

	/// <summary>
	/// 1. Find and store the total mass.
	/// 2. Precompute all accelerations before the first iteration.
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
	int n() const { return (int)tab.size(); }

	/// <summary>
	/// Recall the total mass of particles.
	/// </summary>
	/// <returns></returns>
	double mass() const { return m_mass; }

private:
	/// <summary>
	/// A copy of `tab` made at different times by different
	/// method calls. Only valid within that method call who made it.
	/// </summary>
	V copy;

	/// <summary>
	/// Sum of the masses of all particles.
	/// </summary>
	double m_mass{ 0 };

	/// <summary>
	/// Compute the acceleration felt by particle at index `i` if it were at
	/// the hypothetical location z.
	/// </summary>
	/// <param name="i">Valid index of the particle</param>
	/// <param name="z">The particle's hypothetical location</param>
	/// <returns>Acceleration, or force divided by the particle's mass</returns>
	C accelerate(int i, C const& z) const;
};

