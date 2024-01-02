#include <raylib.hpp>
#include <complex>
#include <deque>

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

/// <summary>
/// A rectangle view-box.
/// </summary>
struct View
{
	/// <summary>
	/// Lesser-coordinate corner (if width and height are positive).
	/// 
	/// In view coordinates.
	/// </summary>
	Cf at;
	/// <summary>
	/// Width and height (normally positive).
	/// 
	/// In view coordinates.
	/// </summary>
	Cf box;
	/// <summary>
	/// Convert from world coordinates to view coordinates
	/// using the axis-aligned scale factor and translation.
	/// </summary>
	/// <param name="w"></param>
	/// <returns></returns>
	Cf w2v(Cf w) const
	{
		float x = w.real() * box.real() + at.real(),
			y = w.imag() * box.imag() + at.imag();
		return Cf(x, y);
	}
	/// <summary>
	/// Convert from view coordinates to world coordinates
	/// using the axis-aligned scale factor and translation.
	/// </summary>
	/// <param name="v"></param>
	/// <returns></returns>
	Cf v2w(Cf v) const
	{
		float x = (v.real() - at.real()) / box.real(),
			y = (v.imag() - at.imag()) / box.imag();
		return Cf(x, y);
	}
};

static void plot(Cf v)
{
	Vector2 v2{ v.real(),v.imag() };
	DrawCircleV(v2, 5.0f, BLACK);
}

int wWinMain(void* _0, void* _1, wchar_t const* _2, int _3)
{
	InitWindow(600, 600, "q2");
	SetTargetFPS(60);

	View view;
	view.box = Cf(300.0f, 300.0f);

	// To use Halton low-discrepancy sequences to fill up the unit square
	// in n-space, use successive prime numbers (2, 3, 5, etc.) and
	// generate x, y, z, etc. coordinates from these individual streams.
	Halton h2(2), h3(3);
	std::deque<Cf> world_halton;

	while (!WindowShouldClose())
	{
		Cf midpoint = Cf(GetScreenWidth() * 0.5f, GetScreenHeight() * 0.5f);
		view.at = midpoint - view.box * 0.5f;

		world_halton.push_front(Cf((float)h2.next(), (float)h3.next()));
		if (world_halton.size() >= 101)
			world_halton.pop_back();

		BeginDrawing();
		{
			ClearBackground(WHITE);

			for (int i = (int)world_halton.size() - 1; i >= 0; i--)
			{
				auto const& world_point = world_halton[i];
				auto const view_point = view.w2v(world_point);
				plot(view_point);
			}

			DrawFPS(16, 16);
		}
		EndDrawing();
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
