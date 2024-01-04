#include "Include.h"

#include <random>

#include "Dyn.h"
#include "Geo2.h"

static Cf c32(C c64) { return Cf((float)c64.real(), (float)c64.imag()); }
static Vector2 v32(Cf c32) { return Vector2{ c32.real(), c32.imag() }; }
static Vector2 v32(C c64) { return v32(c32(c64)); }

static void draw_particle(Dyn const& dyn, int i)
{
	auto color = BLACK;
	auto const& e = dyn[i];
	double prop = e.m / dyn.mass();
	double odds = prop / (1 - prop);
	using std::max;
	using std::min;
	typedef unsigned char U;
	color.a = (U)min(256., odds * 256);
	color.a = max(color.a, (U)50);
#define F(func) func(v32(e.z), (float)e.r, color);
	F(DrawCircleLinesV);
	F(DrawCircleV);
#undef F
}

/// <summary>
/// Universal gravitational constant (units: LLL/T/T/M).
/// </summary>
constexpr double G = 0.05;

/// <summary>
/// Force on the left particle (l) due to the right particle (r).
/// </summary>
/// <param name="l"></param>
/// <param name="r"></param>
/// <returns>Force (units: ML/T/T/T)</returns>
C newton_gravity(Dyn::Entry const& l, Dyn::Entry const& r)
{
	C s = r.z - l.z; double as = abs(s);
	C a = G * r.m * (1 / as / as / as) * s;
	return isfinite(a.real()) && isfinite(a.imag()) ? a : 0;
}

int wWinMain(void* _0, void* _1, void* _2, int _3)
{
	// Simulation (dyn)
	Dyn dyn;
	dyn.par.dt = 0.05;
	{
		// Generate this many (n) particles.
		int constexpr n = 100;
		auto seed = []() { std::random_device dev; return dev(); }();
		auto rng = std::mt19937(seed);
		C const rot = std::polar(1., (double)PI / 6);
		for (int i = n - 1; i >= 0; i--)
		{
			typedef std::uniform_real_distribution<> D;
			D z(-100., 100.), v(-0.1, 0.1), m(1., 5.), r(1., 5.);
			auto sca = [&](D& d) { return d(rng); };
			auto vec = [&](D& d) { return C(sca(d), sca(d)); };
			Dyn::Entry e;
			e.z = vec(z), e.v = vec(v) + rot / abs(e.z) * e.z, e.a = 0;
			e.m = sca(m), e.r = sca(r);
			dyn.tab.push_back(e);
		}
		// It is here where all accelerations are computed
		// for before the first iteration, and where the
		// the total mass (dyn.m_mass) is computed.
		dyn.precompute();
	}
	dyn.drv.pair_force = newton_gravity;

	// Rendering
	float constexpr px_per_l = 1.f;

	// Raylib.
	InitWindow(600, 600, "a");
	SetTargetFPS(60);

	while (!WindowShouldClose())
	{
		dyn.step();
		dyn.bias();

		// The camera allows using the world coordinate system as it is.
		Camera2D cam{};
		cam.offset = Vector2{ GetScreenWidth() / 2.f, GetScreenHeight() / 2.f };
		cam.zoom = px_per_l;

		BeginDrawing();
		{
			ClearBackground(WHITE);
			BeginMode2D(cam);
			for (int i = dyn.n() - 1; i >= 0; i--) draw_particle(dyn, i);
			EndMode2D();
			DrawFPS(16, 16);
		}
		EndDrawing();
	}

	CloseWindow();
	return 0;
}
