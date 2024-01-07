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
		/// Time step (T per step).
		/// </summary>
		double dt{};

		/// <summary>
		/// Inclusive lower and upper bounds for the time step if using variable
		/// time step integrator.
		/// </summary>
		double low_dt{ 0.00005 }, high_dt{ 0.25 };
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
	/// 
	/// No entry is required to exist. If an entry does not exist, some sort of
	/// sensible behavior will be chosen; see the descriptions for individual items.
	/// </summary>
	struct Driver
	{
		/// <summary>
		/// Compute the force on the left (first) particle by the right (second) particle.
		/// 
		/// If this doesn't exist, then the accelerations will either be untouched
		/// after each `step` call or it will be reset to zero.
		/// </summary>
		std::function<C(Entry const&, Entry const&)> pair_force;

		/// <summary>
		/// Judge the strong and weak position vectors calculated in the same run
		/// for some particle in some particle-particle pair to suggest whether to adjust
		/// the time step.
		/// 
		/// Return a positive number if suggested to broaden the time step.
		/// Return negative if suggested to narrow the time step; this answer
		/// also suggests the reaction be retried (best effort) with
		/// a smaller time step. Return 0 if no strong case is made
		/// for either case.
		/// 
		/// If this doesn't exist, then a function that returns "0" for everything
		/// (neutral) is assumed.
		/// </summary>
		std::function<int(C const& strong, C const& weak)> judge_z;

		/// <summary>
		/// See judge_z for information.
		/// 
		/// Summary: Positive number if time step should be increased,
		/// negative if decreased and reaction should be retried as needed, or 0 if
		/// neutral (no suggestion).
		/// 
		/// If not exist: No suggestion will be made.
		/// </summary>
		std::function<int(C const& strong, C const& weak)> judge_v;
	};

	typedef std::vector<Entry> V;

	// For these:
	// Anything except `tab` doesn't make sense to be "moved" without being copied.

	Dyn() = default;
	Dyn(Param const& par) : par(par) {}
	Dyn(Dyn const& dyn) = default;
	Dyn(Dyn&& dyn) noexcept
		: par(dyn.par), tab(std::move(dyn.tab)), drv(dyn.drv), copy()
		, m_mass(dyn.m_mass), m_area(dyn.m_area) {}
	Dyn& operator=(Dyn const& dyn) noexcept
	{
		if (&dyn == this) return *this;
		par = dyn.par, tab = dyn.tab, drv = dyn.drv;
		m_mass = dyn.m_mass, m_area = dyn.m_area;
		copy = V();
		return *this;
	}
	Dyn& operator=(Dyn&& dyn) noexcept
	{
		if (&dyn == this) return *this;
		par = dyn.par, drv = dyn.drv;
		m_mass = dyn.m_mass, m_area = dyn.m_area;
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
	/// 1. Find and store the total mass and area.
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

	/// <summary>
	/// Recall the total "area" of particles.
	/// (All particles have a circular area.)
	/// </summary>
	double area() const { return m_area; }

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
	/// Sum of all areas of all particles.
	/// </summary>
	double m_area{ 0 };

	/// <summary>
	/// Compute the acceleration felt by particle at index `i` if it were at
	/// the hypothetical location z.
	/// </summary>
	/// <param name="i">Valid index of the particle</param>
	/// <returns>Acceleration, or force divided by the particle's mass</returns>
	C accelerate(int i) const;
};

