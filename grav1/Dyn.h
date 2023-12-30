#pragma once

#include "Include.h"
#include <memory>
#include <functional>

/// <summary>
/// An entry in a geometric table.
/// (Stores all particles and their geometric and physical quantities).
/// </summary>
struct DynEntry
{
	/// <summary>
	/// Vectors: position (z: L), velocity (v: L/T), and acceleration (a: L/T/T).
	/// </summary>
	C z, v, a;
	/// <summary>
	/// Scalars: mass (m: M), and the Lennard-Jones "sigma" (ljsigma: L).
	/// </summary>
	double m{}, ljsigma{ 1.0 };
};

/// <summary>
/// Parameters for the simulation.
/// </summary>
struct Param
{
	/// <summary>
	/// The universal gravitational constant (units: LLL / T / T / M)
	/// </summary>
	double g{ 9.0 };
	/// <summary>
	/// Step size (units: T / frame).
	/// </summary>
	double dt{ 0.05 };
	/// <summary>
	/// In the inverse-squared law calculations, a minimum
	/// magnitude of a vector to avoid division by zero or
	/// a number close to zero. (units: L).
	/// </summary>
	double guard0_dist{ 0.30 };
	/// <summary>
	/// An absolute speed limit (units: L/T).
	/// If disabled, INFINITY.
	/// </summary>
	double speed_limit{ 1000000.0 };
	/// <summary>
	/// An absolute acceleration limit (units: L/T/T).
	/// If disabled, INFINITY.
	/// </summary>
	double accel_limit{ 10000000.0 };
	/// <summary>
	/// A constant to be multiplied to the output of
	/// the Lennard-Jones potential gradient
	/// so that when the gradient is multiplied by the distance
	/// the product will have the units of force.
	/// (units: M/T/T).
	/// </summary>
	double lj_force_unit{ 400.0 };
};

/// <summary>
/// A set of behaviors.
/// </summary>
struct DynDriver
{
	/// <summary>
	/// If exists, then a procedure to decide
	/// whether a particle (d) in the table warrants a replacement.
	/// 
	/// Allowed to not exist.
	/// </summary>
	std::function<bool(DynEntry const&)> test_replace{};
	/// <summary>
	/// If exists, then a procedure to generate a new particle
	/// independent of the state of the simulation.
	/// 
	/// Allowed to not exist.
	/// </summary>
	std::function<DynEntry(void)> replace{};
};

class Dyn
{
public:
	int n{};
	Param param;
	DynDriver driver;
	std::shared_ptr<DynEntry[]> tab;
	Dyn(int n, Param param, std::shared_ptr<DynEntry[]> tab) :
		n(n), param(param), tab(tab), mass(compute_mass()) {}
	Dyn(int n, Param param) :
		n(n), param(param), mass(0) {}
	/// <summary>
	/// Compute the acceleration on the particle at this hypothetical location.
	/// </summary>
	C accel(int i, C z) const;
	/// <summary>
	/// Compute the accelerations of all particles and store them.
	/// </summary>
	void accelall() const;
	/// <summary>
	/// Compute the accelerations (see `accel`) and return an abridged copy
	/// of what would be the new row, but don't store it.
	/// </summary>
	DynEntry iter(int i);
	/// <summary>
	/// Advance the state of all particles and store it.
	/// </summary>
	void iterall();
	/// <summary>
	/// Remove the position and velocity bias from all particles.
	/// </summary>
	void center_vz() const;
private:
	/// <summary>
	/// The sum of all masses.
	/// </summary>
	double mass{};
	/// <summary>
	/// Compute the sum of all masses.
	/// </summary>
	/// <returns></returns>
	double compute_mass() const
	{
		double m{};
		for (int i = n - 1; i >= 0; i--) m += tab[i].m;
		return m;
	}
};
