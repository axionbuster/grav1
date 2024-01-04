#include "Include.h"

#include <random>

#include "Dyn.h"
#include "Geo2.h"

static Cf c32(C c64) { return Cf((float)c64.real(), (float)c64.imag()); }
static Vector2 v32(Cf c32) { return Vector2{ c32.real(), c32.imag() }; }
static Vector2 v32(C c64) { return v32(c32(c64)); }

static void draw_particle(Dyn const &dyn, int i)
{
	auto color = BLACK;
	auto const& e = dyn[i];
	color.a = (unsigned char)(e.m / dyn.mass() * 256);
	constexpr unsigned char min_alpha = 50;
	color.a = std::max(color.a, min_alpha);
#define F(func, field) func(v32(e.field), (float)e.r, color);
	F(DrawCircleLinesV, z);
	F(DrawCircleV, z);
#undef F
}

C newton_gravity(Dyn::Param const& par, Dyn::Entry const& l, Dyn::Entry const& r)
{
	C s = r.z - l.z; double as = abs(s);
	C a = par.g * r.m * (1 / as / as / as) * s;
	return isfinite(a.real()) && isfinite(a.imag()) ? a : 0;
}

int wWinMain(void* _0, void* _1, void* _2, int _3)
{
	// Simulation (dyn)
	Dyn dyn;
	dyn.par.g = 3;
	dyn.par.dt = 0.05;
	{
		// Generate this many (n) particles.
		int constexpr n = 50;
		auto seed = []() { std::random_device dev; return dev(); }();
		auto rng = std::mt19937(seed);
		for (int i = n - 1; i >= 0; i--)
		{
			typedef std::uniform_real_distribution<> D;
			D z(-100., 100.), v(-0.1, 0.1), m(1., 5.), r(1., 5.);
			auto sca = [&](D& d) { return d(rng); };
			auto vec = [&](D& d) { return C(sca(d), sca(d)); };
			Dyn::Entry t;
			t.z = vec(z), t.v = vec(v), t.a = 0;
			t.m = sca(m), t.r = sca(r);
			dyn.tab.push_back(t);
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
