#pragma once

#include "Include.h"
#include <memory>

/// <summary>
/// An entry in a kinematical table (a table that
/// stores all kinematic variables but also mass).
/// </summary>
struct KinEntry
{
	/// <summary>
	/// Respectively: position (z), velocity (v), and
	/// acceleration (a); mass (m).
	/// </summary>
	C z, v, a;
	double m{};
};

/// <summary>
/// Parameters for the simulation.
/// </summary>
struct Param
{
	/// <summary>
	/// The universal gravitational constant (units: LLL / T / T / M)
	/// </summary>
	double g{ 10.0 };
	/// <summary>
	/// Step size (units: T).
	/// </summary>
	double dt{ 0.05 };
	/// <summary>
	/// In the inverse-squared law calculations, a minimum
	/// magnitude of a vector to avoid division by zero or
	/// a number close to zero. (units: L).
	/// </summary>
	double guard0_dist{ 0.01 };
	/// <summary>
	/// An absolute speed limit (units: L/T).
	/// </summary>
	double speed_limit{ INFINITY };
};

class Dyn
{
public:
	int n{};
	Param param;
	std::shared_ptr<KinEntry[]> kin;
	Dyn(int n, Param param, std::shared_ptr<KinEntry[]> kin) : n(n), param(param), kin(kin) {}
	Dyn(int n, Param param) : n(n), param(param), kin(new KinEntry[n]) {}
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
	KinEntry iter(int i) const;
	/// <summary>
	/// Advance the state of all particles and store it.
	/// </summary>
	void iterall() const;
	/// <summary>
	/// Compute the barycenter.
	/// </summary>
	/// <returns></returns>
	C bary() const;
	/// <summary>
	/// Remove the position bias from all particles.
	/// </summary>
	void reset_bary() const;
private:
};
