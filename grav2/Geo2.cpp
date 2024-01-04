#include "Geo2.h"

double Halton::next()
{
	x = d - n;
	if (x == 1)
		n = 1, d *= b;
	else
	{
		y = d / b;
		while (x <= y)
			y /= b;
		n = (b + 1) * y - x;
	}
	return n / (double)d;
}

int HaltonQuadrature2D::count(std::function<bool(C const&)> const& test, int n)
{
	int i = n - 1, m = 0;
	while (i-- >= 0) m += !!test(C(h[0].next(), h[1].next()));
	return m;
}

double area_right_crescent(HaltonQuadrature2D& q, C c0, double r0, C c1, double r1, int n)
{
	// Translate both circles so that the left (0) circle is at the origin.
	// Only need to compute that for the right (1) circle, though.
	c1 -= c0;

	// d, m, s, and c are in relative scale (multiples of r0).
	// d --- dimension (side length) of bounding square
	// m --- "midpoint": center of the bounding square
	// s --- "scale": radius of circle 1 (again, as multiples of r0)
	// c --- "center": location of the center of circle 1 (ditto)
	double d{}, m{}, s{}, c{};

	// Let A be c0 as originally provided, and B be c1 (orig.), as well.
	// Let A' be identically the origin, and B' be c1 - c0 (orig.).
	// (Currently, "c1" holds the location of B').
	//
	// In the ray A'B', with A' at the apex, find the "far" (away from A')
	// intersection of the circle A', that point being named P, and
	// the same for the circle B', sim., Q. Make sure the "far" intersection
	// is taken for the point Q.
	//
	// Using the unit circle at the origin A', construct the point 1 + 0i,
	// and call this very point P', because it will be related to Q right now.
	// Note that by a certain choice of a point on the real axis, to be
	// named Q', the two triangles A'P'P and A'Q'Q can be made similar to each other,
	// namely, by making a line parallel to the line segment P'P that
	// passes through Q, and letting the intersection between this line
	// and the real axis be the point Q'.
	//
	// Then, A'P will measure r0 (radius of circle A), A'P' will measure 1,
	// and A'Q will measure the sum (|c1| + r1), and A'Q' will measure
	// (|c1| + r1) / r0. By a similar construction, a line segment of length
	// (|c1| / r0) may be generated as well, likewise on the real axis,
	// and, therefore, the difference r1 / r0 may be extracted.

	// Construct the relative center and radius of the right circle.

	// Construct a circle about A' (origin) passing through B' (true center, circle 2).
	// Then, find the intersection between the circle and the positive real axis (ray).
	// This intersection is given the coordinate (ac1 + 0 i), where i is the imag. unit.
	double ac1 = abs(c1);

	// Now on the ray of the positive real axis, there are four points:
	//  - Radius of circle 0 (r0),
	//  - True center (real) of circle 1 (ac1),
	//  - Radius of circle 1 (r1),
	//  - ac1 extended by the segment (A'r1), distance (ac1 + r1) from A'.
	// Scale everything uniformly so that r0 is mapped to 1, and then
	// by tracing the movements of "ac1" and "r1," find the numbers "s" and "c."
	s = r1 / r0, c = ac1 / r0;

	// Compute the midpoint of the square.
	// First, a bounding rectangle of the two-circle system is created,
	// and its side lengths are measured, and the midpoint (m) taken.
	// The midpoint is on the real axis. By taking the longer edge (length: d)
	// I construct the square.
	//   The left circle is the unit circle at the origin, and the right
	// one is the one we've just computed earlier.
	using std::min;
	using std::max;
	double _l = min(-1., c - s), _r = max(1., c + s);
	double _b = min(-1., -s), _t = max(1., s);
	m = (_l + _r) / 2, d = max(_r - _l, _t - _b);

	// Auxilliary.
	double _s_sq = s * s;

	// The test is conducted by testing each point with the condition:
	// If it belongs to either circle at all,
	// then it belongs to the right circle but not the left. Simple.
	auto const test_inner = [=](C const& p)
		{
			using std::norm; // squared absolute value
			return norm(p) > 1 && norm(p - c) < _s_sq;
		};

	// Since points are given in the unit square at (0, 1) x (0, 1),
	// I will be performing coordinate transformations to fit this
	// to my own boundary circle (centered at (m, 0), with the side length d).
	auto const test = [=](C p)
		{
			p -= C(0.5, 0.5);
			p *= d;
			p += m;
			return test_inner(p);
		};

	// Now, I launch the experiment.
	double f = q.count(test, n) / (double)n;
	double area = (s * s) * (d * d) * f;
	return area;
}
