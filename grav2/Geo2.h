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
	unsigned n{ 0 }, d{ 1 }, x{}, y{}, b{};
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
		: lrr(r0* r0), rrr(r1* r1)
	{
		using std::min;
		using std::max;
		c1 -= c0;
		c = abs(c1);
		double _l = min(-1., c - r1), _r = max(1., c + r1);
		double _b = min(-1., -r1), _t = max(1., r1);
		m = (_l + _r) / 2;
		d = max(_r - _l, _t - _b);
	}

	/// <summary>
	/// Decide whether the point (reoriented coord.) is in the left circle.
	/// </summary>
	bool left(C const& p) const { return std::norm(p) < lrr; }

	/// <summary>
	/// Decide whether the point (reoriented coord.) is in the right circle.
	/// </summary>
	bool right(C const& p) const { return std::norm(p - c) < rrr; }

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
	C from01(C const& h) const { return (h - C(.5, .5)) * d + m; }

	/// <summary>
	/// Recall the area of the bounding square.
	/// </summary>
	double bounding_area() const { return d * d; }

	/// <summary>
	/// Recall the x-coordinate of the center of the bounding square.
	/// </summary>
	double bounding_midpoint_x() const { return m; }
private:
	/// <summary>
	/// Center of the right circle;
	/// Dimension (side length) of the bounding square;
	/// Midpoint of the bounding square;
	/// Squared radius of the left circle (and right circle).
	/// </summary>
	double c{}, d{}, m{}, lrr{}, rrr{};
};
