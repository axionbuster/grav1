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
	float w2v = 75.f;

	InitWindow(600, 600, "Quasi-Monte Carlo Quadrature (WIP)");
	SetTargetFPS(60);

	// Parameters (starting)
	C const orig_c0(0, 1), orig_c1(1, 0);
	double const r0(1.0), r1(2.0);
	int const cap = 100;

	// Angular velocity
	double const omega = 0.005; // radians per frame

	// Frame number
	unsigned fr = 0;

	while (!WindowShouldClose())
	{
		Cf center_px = Cf(GetScreenWidth() * 0.5f, GetScreenHeight() * 0.5f);

		C const c0 = orig_c0 * std::polar(1., omega * fr);
		C const c1 = orig_c1 * std::polar(1., omega * fr);
		GenericLune calculation(c0, r0, c1, r1, cap);

		for (int i = 0; i < cap; i++)
			calculation.lune.advance();

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
						base -= b, tip -= b; // layout bug
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
				auto const co = calculation.lune.in(p) ? DARKGREEN : BLACK;
				plot(at, co);
			}

			// Show statistics.
			DrawFPS(16, 16);

			fr++;
		}
		EndDrawing();
	}

	CloseWindow();
	return 0;
}
