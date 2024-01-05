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
	double prop = e.m / dyn.mass(); // kinda doesn't make sense...
	using std::max;
	using std::min;
	typedef unsigned char U;
	color.a = (U)min(256., prop * 256);
	color.a = max(color.a, (U)50);
#define F(func) func(v32(e.z), (float)e.r, color);
	F(DrawCircleLinesV);
	F(DrawCircleV);
#undef F
}

/// <summary>
/// Universal gravitational constant (units: LLL/T/T/M).
/// </summary>
constexpr double G = 1.;

/// <summary>
/// Decide whether both components of the complex number are finite.
/// </summary>
static bool finite(C const& c) { return isfinite(c.real()) && isfinite(c.imag()); }

/// <summary>
/// Force on the left particle (l) due to the right particle (r).
/// </summary>
/// <param name="l"></param>
/// <param name="r"></param>
/// <returns>Force (units: ML/T/T/T)</returns>
C newton_gravity(Dyn::Entry const& l, Dyn::Entry const& r)
{
	// Low-discrepancy sequences needed for parts of the calculation.
	static Halton2D hh;

	C s = r.z - l.z; double as = abs(s);

	if (as < l.r + r.r)
	{
		// The circles representing them intersect.
		// The simple calculation below doesn't apply.
		// So, integrate the infinitesimal forces to get the total force for each
		// small patch of the region of the left circle that is outside
		// the right circle.

		// Number of Monte-Carlo trials.
		constexpr int M = 100;
		// Force per mass (integrated).
		C fpm;
		// Number of pieces sampled in the left crescent (one-sided lune---just "lune").
		int n{};
		// Circular intersection.
		CircularIntersection sect(l.z, l.r, r.z, r.r);
		// Computation of lunar area and force.
		auto compute = [&](C const& p)
			{
				// `p` is sampled from a certain square in
				// a reoriented coordinate system, where
				// the left circle is centered at the origin, and
				// the right circle as at (`as`, 0). Lengths have
				// not changed.

				if (!sect.left(p) || sect.right(p)) return;
				n++;
				C arm = as - p; double dist = abs(arm);
				// [***] `fpm` will be missing the factors of: M (=# trials), G, dm.
				fpm += 1 / dist / dist / dist / M * sect.orient(arm);
			};
		// Boom.
		for (int i = M - 1; i >= 0; i--)
			sect.monte(hh.next(), compute);
		// If no sample hit, the integration failed.
		if (!n) return 0;
		// Compute dm by the ratio -> dm : m = 1 : n
		// (where dm: infinitesimal mass, m: mass of left particle,
		// n: number of samples hit).
		double dm = l.m / n;
		// [***] Multiply back the missing factors in `fpm`.
		fpm *= M * G * dm;
		return finite(fpm) ? fpm : 0;
	}
	else
	{
		C f = G * l.m * r.m * (1 / as / as / as) * s;
		return finite(f) ? f : 0;
	}
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
		C const rot = std::polar(0., PI64 / 6);
		for (int i = n - 1; i >= 0; i--)
		{
			std::uniform_real_distribution<> z(-50, 50.), v(-1, 1);
			std::cauchy_distribution<> m(5., 1.), r(5., 1.); // center; scale.
#define sca(d) d(rng)
#define vec(d) C(sca(d), sca(d))
			Dyn::Entry e;
			e.z = vec(z), e.v = vec(v) + rot / abs(e.z) * e.z, e.a = 0;
			e.m = abs(sca(m)), e.r = abs(sca(r));
#undef vec
#undef sca
			dyn.tab.push_back(e);
		}
		// It is here where all accelerations are computed
		// for before the first iteration, and where the
		// the total mass (dyn.m_mass) is computed.
		dyn.precompute();
	}
	dyn.drv.pair_force = newton_gravity;

	// Rendering
	float constexpr px_per_l = 2.f;

	// Raylib.
	InitWindow(600, 600, "Gravity");
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
