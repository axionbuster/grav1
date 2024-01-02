#include <raylib.hpp>
#include <complex>
#include <deque>

using namespace std::literals::complex_literals;

typedef std::complex<double> C;
typedef std::complex<float> Cf;

/// <summary>
/// Generate a Halton sequence (algorithm is due to Wikipedia) of a given base (`b`).
/// </summary>
struct Halton
{
	int n{ 0 }, d{ 1 }, x{}, y{}, b{};
	Halton(int b) : b(b) {}
	double next();
};

static void plot(Cf v)
{
	Cf const dim(5.f, 5.f);
	v -= dim * 0.5f;
	Vector2 v2{ v.real(), v.imag() },
		dim2{ dim.real(), dim.imag() };
	DrawRectangleV(v2, dim2, BLACK);
}

static void line(Cf from, Cf to, float thick = 1.f, Color c = BLACK)
{
	Vector2 from2{ from.real(), from.imag() },
		to2{ to.real(), to.imag() };
	DrawLineEx(from2, to2, thick, c);
}

static void axis1(Cf v, Cf interval, int n, char const* format, double show_unit = 0.f)
{
	Cf const tick = 5.if;
	Cf const rot = interval / abs(interval);
	Cf const neg_extr = -interval * (float)n,
		pos_extr = interval * (float)n;
	line(v + neg_extr, v + pos_extr);

	for (float i = (float)-n; i <= n; i++)
	{
		if (i == 0) continue;
		Cf const at = interval * i,
			to1 = at + rot * tick * .5f,
			to2 = at - rot * tick * .5f;
		line(v + to1, v + to2);

		if (show_unit == 0) continue;
		// Some font stuff...
		Font font = GetFontDefault();
		float font_size =
			(i == 1 || i == -1) ? 16.f : 10.f;
		float spacing = 1.0f;
		Color tint = BLACK;

		// Plan printing.
		char text[256];
		snprintf(text, sizeof(text), format, i * show_unit);
		Vector2 text_dim = MeasureTextEx(font, text, font_size, spacing);
		Cf text_dim_cf(text_dim.x, text_dim.y);

		// Center.
		Cf text_loc = at - text_dim_cf * .5f,
			text_loc_v = v + text_loc;
		Vector2 text_loc_v2{ text_loc_v.real(), text_loc_v.imag() };

		// Print.
		DrawTextEx(font, text, text_loc_v2, font_size, spacing, tint);
	}
}

static void axes(Cf v, Cf re, Cf im, int n, double label = 0)
{
	axis1(v, re, n, "%.2f", label);
	axis1(v, im, n, "%.2fi", label);
}

static void circle(Cf o, float r, Color const* fill, Color const* border)
{
	Vector2 v{ o.real(), o.imag() };
	if (fill)
		DrawCircleV(v, r, *fill);
	if (border)
		DrawCircleLinesV(v, r, *border);
}

int wWinMain(void* _0, void* _1, wchar_t const* _2, int _3)
{
	// Dimensions of the unit square in world coordinates in view coordinates.
	// (pixels per L).
	float w2v = 100.f;

	InitWindow(600, 600, "q2");
	SetTargetFPS(60);

	// To use Halton low-discrepancy sequences to fill up the unit square
	// in n-space, use successive prime numbers (2, 3, 5, etc.) and
	// generate x, y, z, etc. coordinates from these individual streams.
	Halton h2(2), h3(3);
	std::deque<Cf> world_halton;

	// Test rotating axes.
	unsigned fr = 0;

	while (!WindowShouldClose())
	{
		// Compute the midpoint.
		Cf midpoint = Cf(GetScreenWidth() * 0.5f, GetScreenHeight() * 0.5f);

		// Update the low-discrepancy square.
		world_halton.push_front(Cf((float)h2.next(), (float)h3.next()));
		if (world_halton.size() >= 101)
			world_halton.pop_back();

		BeginDrawing();
		{
			ClearBackground(WHITE);

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

			// Plot the low-discrepancy square.
			for (int i = (int)world_halton.size() - 1; i >= 0; i--)
			{
				auto const& world_point = world_halton[i];
				auto const view_point = midpoint + Cf(world_point.real() * w2v,
					world_point.imag() * w2v);
				plot(view_point);
			}

			DrawFPS(16, 16);
		}
		EndDrawing();

		fr++;
	}

	CloseWindow();
	return 0;
}

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
