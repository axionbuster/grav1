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

struct Lune
{
	Halton h2, h3;
	std::deque<C> log;
	double c{}, rsq{};
	int freq{}, cap{};
	Lune(Halton h2, Halton h3, double c, double r, int cap = 100)
		: h2(h2), h3(h3), c(c), rsq(r* r), cap(cap)
	{
		using std::min;
		using std::max;
		// Starting from a rectangle, construct the bounding square.
		double le = min(-1., c - r),
			ri = max(1., c + r),
			to = max(1., r),
			bo = min(-1., -r);
		dim = max(max(abs(le), abs(ri)),
			max(abs(to), abs(bo))) * 2;
	}
	void advance();
	double quadrature() const
	{
		return dim * dim * freq / log.size();
	}
	double dimension() const { return dim; }
	static bool left_static(C p) { return std::norm(p) < 1.0; }
	bool left(C p) const { return left_static(p); }
	bool right(C p) const { return std::norm(p - c) < rsq; }
	bool in(C p) const { return left(p) && !right(p); }
private:
	double dim{};
};

int wWinMain(void* _0, void* _1, wchar_t const* _2, int _3)
{
	// Dimensions of the unit square in world coordinates in view coordinates.
	// (pixels per L).
	float w2v = 100.f;

	InitWindow(600, 600, "q2");
	SetTargetFPS(60);

	// As a test case, use the Lune of Hippocrates.
	Lune lune(
		// Use successive prime numbers as the bases for the Halton sequences.
		Halton(2), Halton(3),
		// Center of the other circle; radius thereof.
		1 / sqrt(2), 1 / sqrt(2),
		// Capacity (number of evaluations to hold in queue).
		1000
	);

	while (!WindowShouldClose())
	{
		lune.advance();
		Cf midpoint = Cf(GetScreenWidth() * 0.5f, GetScreenHeight() * 0.5f);

		BeginDrawing();
		{
			ClearBackground(WHITE);

			// Some decor...
			circ_axes(w2v, midpoint);

			// The chosen square...
			{
				Rectangle r{};
				float w2vf = w2v;
				r.x = (float)lune.dimension() / -2 * w2vf + (float)midpoint.real();
				r.y = (float)lune.dimension() / -2 * w2vf + (float)midpoint.imag();
				r.width = (float)lune.dimension() * w2vf;
				r.height = (float)lune.dimension() * w2vf;
				DrawRectangleLinesEx(r, 1.0f, BLACK);
			}

			// Plotting.
			for (int i = (int)(lune.log.size() - 1); i >= 0; i--)
			{
				auto const& p = lune.log[i];
				auto const at = downgrade(p) * w2v + midpoint;
				Color c = BLACK;
				if (lune.in(p))
					c = DARKGREEN;
				plot(at, c);
			}

			// Statistics.
			DrawFPS(16, 16);

			char msg[512];
			snprintf(msg, sizeof msg, "rel freq: %.4f; quadrature: %.4f\n"
				"size: %d; cap: %d.",
				(double)lune.freq / lune.log.size(), lune.quadrature(),
				(int)lune.log.size(), (int)lune.cap
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
	c1fill.a = 100;
	Color c2fill = BLUE,
		c2border = BLUE;
	c2fill.a = 100;
	float const circle2_radius = 1 / sqrtf(2);
	circle(midpoint, 1.f * w2v, &c1fill, &c1border);
	circle(midpoint + w2v * circle2_radius, w2v * circle2_radius, &c2fill, &c2border);
}

void Lune::advance()
{
	// Construct a point in the unit square in (0,1) x (0,1),
		// then scale and translate it into the bounding square.
	C p0 = C(h2.next(), h3.next()) - C(0.5, 0.5),
		p = p0 * dim;

	log.push_back(p);
	if (log.size() > cap)
	{
		if (in(log.front()))
			freq--;
		log.pop_front();
	}

	if (in(p))
		freq++;
}
