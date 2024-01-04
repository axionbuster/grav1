// Demo. (Under major restructuring right now.)

#include <deque>

#include "Header.h"
#include "Q2vis.h"
#include "Lune.h"

using namespace vis;
using namespace lune;

/// <summary>
/// The entry point (Windows).
/// </summary>
int wWinMain(void* _0, void* _1, wchar_t const* _2, int _3)
{
	// Dimensions of the unit square in world coordinates in view coordinates.
	// (pixels per L).
	float w2v = 25.f;

	InitWindow(600, 600, "Quasi-Monte Carlo Quadrature (WIP)");
	SetTargetFPS(60);

	// Parameters (starting)
	C const orig_c0(0, 1), orig_c1(1, 0);
	double const r0(1.0), r1(2.0);
	int const cap = 100;

	// Angular velocity
	double const spin = 0.005; // radians per frame
	double const orbit = 0.0025; // radians per frame
	double const orbit_arm_w = 4; // in units of world coordinates

	// Frame number
	unsigned fr = 0;

	// For the sources of "randomness" (in fact not random for fast convergence),
	// pick up where I left off.
	halton::Halton h2(2), h3(3);

	// Smoothing for statistical reporting.
	struct Stat
	{
		double relfreq{}, quadrature{};
	};
	std::deque<Stat> stats;
	int const stats_cap = 500;
	auto const put_stat = [&](Stat const& s)
		{
			if (stats.size() >= stats_cap) stats.pop_front();
			stats.push_back(s);
		};
	struct StatSummary
	{
		double mean_relfreq{}, s_stdev_relfreq{};
		double mean_quadrature{}, s_stdev_quadrature{};
	};
	auto const summarize = [&]() -> StatSummary
		{
			// Apply Welford's online algorithm for the calculation
			// of summary statistics.
			int i;
			double m1r{}, m2r{}, m1q{}, m2q{}; // 1st and 2nd moments
			for (i = 1; i <= stats.size(); i++)
			{
				double r = stats[i - 1].relfreq, q = stats[i - 1].quadrature;
				double m1r0 = m1r, m1q0 = m1q;
				m1r += (r - m1r) / i;
				m1q += (q - m1q) / i;
				m2r += (r - m1r0) * (r - m1r);
				m2q += (q - m1q0) * (q - m1q);
			}
			StatSummary s{};
			s.mean_relfreq = m1r, s.mean_quadrature = m1q;
			if (i) s.s_stdev_relfreq = m2r / (i - 1), s.s_stdev_quadrature = m2q / (i - 1);
			return s;
		};

	while (!WindowShouldClose())
	{
		using std::move;

		Cf center_px = Cf(GetScreenWidth() * 0.5f, GetScreenHeight() * 0.5f);

		C const c0 = orig_c0 * std::polar(1., spin * fr) + std::polar(orbit_arm_w, orbit * fr);
		C const c1 = orig_c1 * std::polar(1., spin * fr) + std::polar(orbit_arm_w, orbit * fr);
		GenericLune calculation(c0, r0, c1, r1, cap);

		calculation.lune.swap_halton(move(h2), move(h3));
		for (int i = 0; i < cap; i++)
			calculation.lune.advance();
		calculation.lune.swap_halton(move(h2), move(h3));

		BeginDrawing();
		{
			ClearBackground(WHITE);

			// Crosshair at the origin.
			{
				Cf a = 5, b = 1.f;
				for (int i = 0; i < 4; i++)
				{
					Cf base = center_px + b, tip = base + a;
					if (i == 0 || i == 1)
						base -= b, tip -= b; // deal with layout bug of Raylib
					DrawLineV(v2(base), v2(tip), BLACK);
					a *= 1.if, b *= 1.if;
				}
			}

			// Draw the bounding rectangle.
			{
				auto rec = calculation.bounding();
				rec.transform(w2v, center_px);
#define side(m, n) DrawLineV(v2(rec.c[m]), v2(rec.c[n]), BLACK)
				side(0, 1), side(1, 2), side(2, 3), side(3, 0);
#undef side
			}

			// The circles.
			{
				Cf c0_prime = downgrade(c0) * w2v + center_px;
				Cf c1_prime = downgrade(c1) * w2v + center_px;
				float r0_prime = (float)r0 * w2v, r1_prime = (float)r1 * w2v;
				auto color0 = RED, color1 = BLUE;
				constexpr auto alpha256 = 100;
				color0.a = alpha256, color1.a = alpha256;
#define circle(n) DrawCircleLinesV(v2(c##n##_prime), r##n##_prime, color##n), \
DrawCircleV(v2(c##n##_prime), r##n##_prime, color##n)
				circle(0), circle(1);
#undef circle
			}

			// The points.
			for (int i = (int)(calculation.lune.log.size() - 1); i >= 0; i--)
			{
				auto const& p = calculation.lune.log[i];
				auto const at = downgrade(calculation.invert(p)) * w2v + center_px;
				auto const co = calculation.lune.in(p) ? YELLOW : BLACK;
				plot(at, co);
			}

			// Compute and show statistics.
			{
				Stat s{};
				s.relfreq = (double)calculation.lune.freq / calculation.lune.log.size();
				s.quadrature = abs(calculation.homt) * calculation.lune.quadrature();
				put_stat(s); // Influences following `summarize` call.
			}
			StatSummary summary = summarize();

			DrawFPS(16, 16);
			char msg[999];
			snprintf(msg, sizeof(msg), "%d-frame statistics\n"
				"relfreq\n\tmean: %.3f\n\tstdev: %.5f\n"
				"quadrature\n\tmean: %.3f\n\tstdev: %.3f\n"
				"(sample stdev; each frame)",
				(int)stats.size(),
				summary.mean_relfreq, summary.s_stdev_relfreq,
				summary.mean_quadrature, summary.s_stdev_quadrature);
			DrawText(msg, 16, 48, 20, BLACK);

			fr++;
		}
		EndDrawing();
	}

	CloseWindow();
	return 0;
}
