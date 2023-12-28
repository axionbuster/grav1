#pragma once

#include "Include.h"
#include <memory>

using namespace std::literals::complex_literals;

struct KinEntry
{
	C z, v, a;
	double m{};
};

struct Param
{
	double g{ 0.01 }, dt{ 0.05 }, guard0_dist{ 0.0001 }, speed_limit{ 100.0 };
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
private:
};
