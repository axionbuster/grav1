#include <deque>

#include "Header.h"
#include "Q2vis.h"
#include "Halton.h"

using namespace vis;
using namespace halton;

/// <summary>
/// Plot the circles and the axes.
/// </summary>
/// <param name="w2v"></param>
/// <param name="midpoint"></param>
void circ_axes(float w2v, Cf midpoint);

/// <summary>
/// Interactive quadrature of the area of the one-sided lune (crescent).
/// 
/// When a circle with a given center on the x-axis and radius
/// intersects with the unit circle at the origin, find the area
/// of the circle on the right except that part which belongs to
/// the central unit circle, if that intersecting region exists (or else, zero).
/// </summary>
struct Lune
{
	/// <summary>
	/// Sequence of points that have been recently sampled (FIFO).
	/// 
	/// Back = recent, front = later.
	/// </summary>
	std::deque<C> log;
	/// <summary>
	/// x-coordinate of the center of the right circle;
	/// squared radius of the right circle.
	/// 
	/// The left circle is always the unit circle (radius 1) at the origin.
	/// </summary>
	double c{}, rsq{};
	/// <summary>
	/// The number of points satisfying the criteria;
	/// number of points to keep in the log.
	/// </summary>
	int freq{}, cap{};
	/// <summary>
	/// Set up the computation of the area of the right-side of the lune
	/// (crescent) by the intersection of the unit circle at the origin
	/// and a second circle centered at (c, 0) with r > 0 as the radius.
	/// </summary>
	/// <param name="c">x-coordinate of the center of the second circle.</param>
	/// <param name="r">Positive radius of the second circle.</param>
	/// <param name="cap">Number of points to hold in history queue (`log`).</param>
	Lune(double c, double r, int cap = 100)
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
	/// <summary>
	/// Sample a point and update internal statistics.
	/// </summary>
	void advance();
	/// <summary>
	/// Once at least one point has been sampled,
	/// compute the quadrature by inspection of the internal statistics.
	/// 
	/// (Simple arithmetic).
	/// </summary>
	/// <returns></returns>
	double quadrature() const
	{
		return dim * dim * freq / log.size();
	}
	/// <summary>
	/// Recall the side length of the bounding square.
	/// </summary>
	/// <returns></returns>
	double dimension() const { return dim; }
	/// <summary>
	/// Recall the x-coordinate of the center of the bounding square.
	/// </summary>
	/// <returns></returns>
	double xmidpoint() const { return m_xmidpoint; }
	/// <summary>
	/// Decide whether the point belongs to the left circle.
	/// </summary>
	/// <param name="p"></param>
	/// <returns></returns>
	static bool left_static(C p) { return std::norm(p) < 1.0; }
	/// <summary>
	/// Decide whether the point belongs to the left circle.
	/// </summary>
	/// <param name="p"></param>
	/// <returns></returns>
	bool left(C p) const { return left_static(p); }
	/// <summary>
	/// Decide whether the point belongs to the right circle.
	/// </summary>
	/// <param name="p"></param>
	/// <returns></returns>
	bool right(C p) const { return std::norm(p - c) < rsq; }
	/// <summary>
	/// Decide whether the point belongs to the right circle
	/// but not the left circle (used for quadrature).
	/// </summary>
	/// <param name="p"></param>
	/// <returns></returns>
	bool in(C p) const { return !left(p) && right(p); }
private:
	/// <summary>
	/// Internal low-discrepancy sequences used to evenly generate
	/// points in the constructed bounding square, of successive
	/// prime numbers as the "bases" (peculiar details of the algorithm).
	/// 
	/// For each "random" point generated with these:
	/// First a term at h2 is chosen (0 &lt; [h2] &lt; 1), and is turned into the
	/// x-coordinate of that point. The same is done for h3, which becomes
	/// the y-coordinate of that point.
	/// </summary>
	Halton h2, h3;
	/// <summary>
	/// The positive side length of the bounding square.
	/// </summary>
	double dim{};
	/// <summary>
	/// Center of the bounding square.
	/// </summary>
	double m_xmidpoint{};
};

/// <summary>
/// The entry point (Windows).
/// </summary>
int wWinMain(void* _0, void* _1, wchar_t const* _2, int _3)
{
	// Dimensions of the unit square in world coordinates in view coordinates.
	// (pixels per L).
	float w2v = 100.f;

	InitWindow(600, 600, "q2");
	SetTargetFPS(60);

	// As a test case, use the Lune of Hippocrates.
	Lune hippocrates(
		// Center of the other circle; radius thereof.
		1 / sqrt(2), 1 / sqrt(2),
		// Capacity (number of evaluations to hold in queue).
		1000
	);

	while (!WindowShouldClose())
	{
		// Do a few times.
		for (int i = 0; i < 100; i++)
			hippocrates.advance();
		Cf midpoint = Cf(GetScreenWidth() * 0.5f, GetScreenHeight() * 0.5f);

		BeginDrawing();
		{
			ClearBackground(WHITE);

			// Some decor...
			circ_axes(w2v, midpoint);

			// The chosen square (outlines)...
			{
				Rectangle r{};
				float w2vf = w2v;
				r.x = (float)(hippocrates.dimension() / -2
					+ hippocrates.xmidpoint()) * w2vf + (float)midpoint.real();
				r.y = (float)hippocrates.dimension() / -2 * w2vf + (float)midpoint.imag();
				r.width = (float)hippocrates.dimension() * w2vf;
				r.height = (float)hippocrates.dimension() * w2vf;
				DrawRectangleLinesEx(r, 1.0f, BLACK);
			}

			// Plotting.
			for (int i = (int)(hippocrates.log.size() - 1); i >= 0; i--)
			{
				auto const& p = hippocrates.log[i];
				auto const at = downgrade(p) * w2v + midpoint;
				Color c = BLACK;
				if (hippocrates.in(p))
					c = DARKGREEN;
				plot(at, c);
			}

			// Statistics.
			DrawFPS(16, 16);
			char msg[512];
			snprintf(msg, sizeof msg, "rel freq: %.4f\nquadrature: %.4f\n"
				"size: %d\ncap: %d\n"
				"dim: %.4f",
				(double)hippocrates.freq / hippocrates.log.size(), hippocrates.quadrature(),
				(int)hippocrates.log.size(), (int)hippocrates.cap,
				hippocrates.dimension()
			);
			DrawText(msg, 16, 40, 20, BLACK);
		}
		EndDrawing();
	}

	CloseWindow();
	return 0;
}

void circ_axes(float w2v, Cf midpoint)
{
	// Plot the axes.
	axes
	(
		midpoint,
		w2v,
		w2v * 1.if,
		// How many tick marks?
		1 + (int)(midpoint.real() / w2v),
		// Scale factor.
		1.f
	);

	// Plot the circles.
	// (Setting up the lunes of Hippocrates).
	Color c1fill = RED,
		c1border = RED;
	c1fill.a = 100; // Transparency (low = transparent; proportion out of 256 parts).
	Color c2fill = BLUE,
		c2border = BLUE;
	c2fill.a = 100;
	// Radius and position of the second circle is hard-coded for the case of
	// the lune of Hippocrates.
	float const circle2_radius = 1 / sqrtf(2);
	circle(midpoint, 1.f * w2v, &c1fill, &c1border);
	circle(midpoint + w2v * circle2_radius, w2v * circle2_radius, &c2fill, &c2border);
}

void Lune::advance()
{
	// Construct a point in the unit square in (0,1) x (0,1),
	// then scale and translate it into the bounding square.
	C p0 = C(h2.next(), h3.next()) - C(0.5, 0.5),
		p = p0 * dim + m_xmidpoint;

	// Logging, not important for the computation.
	log.push_back(p);
	if (log.size() > cap)
	{
		if (in(log.front()))
			freq--;
		log.pop_front();
	}

	// Monte Carlo.
	if (in(p))
		freq++;
}
