#pragma once

#include "Include.h"

/// <summary>
/// Generate a Halton sequence (algorithm is due to Wikipedia) of a given base (`b`).
/// 
/// To use Halton low-discrepancy sequences to fill up the unit square
/// in n-space, use successive prime numbers (2, 3, 5, etc.) and
/// generate x, y, z, etc. coordinates from these individual streams.
/// </summary>
struct Halton
{
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
private:
	unsigned n{ 0 }, d{ 1 }, x{}, y{}, b{};
	friend struct Halton2D;
	/// <summary>
	/// Create a default instance in an unspecified (erroneous) state.
	/// </summary>
	Halton() : b(0) {}
};

/// <summary>
/// Using the Halton low-discrepancy sequences,
/// generates points evenly in the unit square in (0,1) x (0,1).
/// </summary>
struct Halton2D
{
	/// <summary>
	/// Exposed internal state (Halton low-discrepancy sequences).
	/// </summary>
	Halton h[2];

	/// <summary>
	/// Create a valid instance of Halton2D.
	/// </summary>
	Halton2D() { h[0] = Halton(2), h[1] = Halton(3); }

	/// <summary>
	/// Generate a point in the (0,1) x (0,1) square.
	/// </summary>
	C next() { return C(h[0].next(), h[1].next()); }
};

/// <summary>
/// Store a pair of possibly intersecting circles in a reoriented coordinate system,
/// where the left circle is centered at the origin, and the right circle
/// has a center on the positive real axis. Here lengths do not change;
/// only orientation does.
/// </summary>
struct CircularIntersection
{
public:
	/// <summary>
	/// Construct data about the (possible) intersection between two circles.
	/// </summary>
	/// <param name="c0">Center of the left circle.</param>
	/// <param name="r0">Radius of the left circle.</param>
	/// <param name="c1">Center of the right circle.</param>
	/// <param name="r1">Radius of the right circle.</param>
	CircularIntersection(C c0, double r0, C c1, double r1)
		: lr(r0), lrsq(r0* r0), rrsq(r1* r1)
	{
		// Translate (geometry) as required.
		c1 -= c0;
		// Make a circle centered about (0,0) passing through point c1
		// and then construct the intersection (c) between this circle
		// and the positive real axis ray. This represents rotation.
		c = abs(c1);
		// Compute the reverse rotation.
		derot = c1 / c;
	}

	/// <summary>
	/// Decide whether the point (reoriented coord.) is in the left circle.
	/// </summary>
	bool left(C const& p) const { return std::norm(p) < lrsq; }

	/// <summary>
	/// Decide whether the point (reoriented coord.) is in the right circle.
	/// </summary>
	bool right(C const& p) const { return std::norm(p - c) < rrsq; }

	/// <summary>
	/// Perform a trial, calling the user-defined process `f` that takes in
	/// the point `p` which is in the reoriented coordinate system.
	/// </summary>
	/// <param name="h">A point in the (0,1) x (0,1) square</param>
	/// <param name="f">A process for which p is in the reoriented coordinate system</param>
	void monte(C const& h, std::function<void(C const& p)> const& f) const { f(from01(h)); }

	/// <summary>
	/// Transform a point in the (0,1) x (0,1) square to the bounding square.
	/// </summary>
	C from01(C const& h) const { return (h - C(.5, .5)) * lr; }

	/// <summary>
	/// Orient the reoriented vector to the original orientation (but the
	/// left circle will be still at the origin). Essentially, un-rotate.
	/// </summary>
	C unrotate(C const& p) const { return p * derot; }
private:
	/// <summary>
	/// Center of the right circle (x-coordinate);
	/// Midpoint of the bounding square;
	/// Radius of the left circle;
	/// Squared radius of the left circle (and right circle).
	/// </summary>
	double c{}, lr{}, lrsq{}, rrsq{};
	/// <summary>
	/// Rotation needed to transform the reoriented coordinate system vector
	/// to the original coordinate system (but without translation).
	/// 
	/// "De-rotation."
	/// </summary>
	C derot;
};
