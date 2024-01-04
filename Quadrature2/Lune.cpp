#include "Lune.h"

using namespace lune;

Lune::Lune(double c, double r, int cap)
	: h2(2), h3(3), c(c), rsq(r* r), cap(cap)
{
	using std::min;
	using std::max;
	// Locate the medians (midpoints) of a bounding rectangle.
	double le = min(-1., c - r),
		ri = max(1., c + r),
		to = max(1., r),
		bo = min(-1., -r);
	// Construct the midpoint of the left and right medians
	// to center the bounding square.
	m_xmidpoint = (le + ri) / 2;
	// Compute the side length of a safe bounding square.
	dim = max(ri - le, to - bo);
}

void Lune::advance()
{
	// Construct a point in the unit square in (0,1) x (0,1),
	// then scale and translate it into the bounding square.
	C p0 = C(h2.next(), h3.next()) - C(0.5, 0.5);
	C p = p0 * dim + m_xmidpoint;

	// Logging, not important for the computation.
	log.push_back(p);
	if (log.size() > cap)
	{
		if (in(log.front())) freq--;
		log.pop_front();
	}

	// Monte Carlo.
	if (in(p))
		freq++;
}

void GenericLune::BoundingSquare::transform(C homt, C tr)
{
#define par(op, x) c[0] op x, c[1] op x, c[2] op x, c[3] op x
	par(*=, homt);
	par(+=, tr);
#undef par
}

GenericLune::GenericLune(C c0, double r0, C c1, double r1, int cap)
	: c0(c0), r0(r0), c1(c1), r1(r1), lune(0, 0)
{
	c1 -= c0;
	double const ac1 = abs(c1);
	tr = c0;
	homt = r0 / ac1 * c1;
	lune = Lune(ac1 / r0, r1 / r0, cap);
}

GenericLune::BoundingSquare GenericLune::bounding() const
{
	BoundingSquare sq
	{
		C(-0.5, -0.5), C(0.5, -0.5),
		C(0.5, 0.5), C(-0.5, 0.5)
	};
	// Homothety (scaling and rotation), and then translation (complex).
	sq.transform(lune.dimension(), lune.xmidpoint());
	sq.transform(homt, tr);
	return sq;
}

C GenericLune::invert(C point) const
{
	return point * homt + tr;
}
