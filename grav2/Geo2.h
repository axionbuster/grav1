#pragma once

#include "Include.h"

#include <functional>

/// <summary>
/// Generate a Halton sequence (algorithm is due to Wikipedia) of a given base (`b`).
/// 
/// To use Halton low-discrepancy sequences to fill up the unit square
/// in n-space, use successive prime numbers (2, 3, 5, etc.) and
/// generate x, y, z, etc. coordinates from these individual streams.
/// </summary>
struct Halton
{
	unsigned n{ 0 }, d{ 1 }, x{}, y{}, b{};
public:
	/// <summary>
	/// Initialize a sequence with the given prime base.
	/// </summary>
	/// <param name="b">Base (prime number).</param>
	Halton(unsigned b) : b(b) {}
	/// <summary>
	/// Extract a term and advance the internal state.
	/// </summary>
	/// <returns>A number in the open interval (0, 1).</returns>
	double next();
};

/// <summary>
/// Quasi-Monte Carlo integrator in two dimensions
/// using Halton low-discrepancy sequences for fast convergence.
/// 
/// Region of integration: Unit square in (0, 1) x (0, 1).
/// </summary>
struct HaltonQuadrature2D
{
	/// <summary>
	/// Exposed internal state (Halton low-discrepancy sequences).
	/// </summary>
	Halton h[2];
	/// <summary>
	/// Compute the frequency of the points that pass the test.
	/// 
	/// Points are generated using the internal array of sequences (`haltons`),
	/// whose states will change.
	/// </summary>
	/// <param name="test">A test against a point in x:(0, 1) by y:(0, 1), in this order,
	/// expressed as a complex number (real is x, imag is y).</param>
	/// <param name="n">Positive sample count.</param>
	/// <returns>Frequency of `true` values returned by the test in the closed interval [0, n].</returns>
	int count(std::function<bool(C const &)> const& test, int n = 100);
};

/// <summary>
/// Compute the area of the crescent in the right disk (circle 1)
/// formed by the intersection between that and the left disk (circle 0),
/// if such an intersection exists; otherwise, 0.
/// </summary>
/// <param name="q">Halton sequence bundle.</param>
/// <param name="c0">Center of the left circle.</param>
/// <param name="r0">Positive radius of the left circle.</param>
/// <param name="c1">Center of the right circle.</param>
/// <param name="r1">Positive radius of the right circle.</param>
/// <param name="n">Positive number of trials for the quasi-Monte Carlo quadrature.</param>
/// <returns>The approximate area, or 0.</returns>
double area_right_crescent(
	HaltonQuadrature2D &q,
	C c0, double r0,
	C c1, double r1,
	int n = 100
);
