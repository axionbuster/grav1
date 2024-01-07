#include "Include.h"

#include <random>

#include "Dyn.h"
#include "Geo2.h"

using namespace dyn;

static Cf c32(C c64) { return Cf((float)c64.real(), (float)c64.imag()); }
static Vector2 v32(Cf c32) { return Vector2{ c32.real(), c32.imag() }; }
static Vector2 v32(C c64) { return v32(c32(c64)); }

static void draw_particle(Dyn const& dyn, int i)
{
	auto color = BLACK;
	auto const& e = dyn[i];
	auto circarea = [](double r) { return r * r * PI64; };
	double score = (e.m / dyn.mass()) / (circarea(e.r) / dyn.area());
	score = score / (1 + score);
	using std::max;
	using std::min;
	typedef unsigned char U;
	color.a = (U)min(250., score * 256);
	color.a = max(color.a, (U)50);
	Vector2 z_postproc;
	float r_postproc;
	{
		// Squish.
		double r = abs(e.z);
		double r2 = 250. * tanh(r / 250.);
		double ratio = r2 / r;
		z_postproc = v32(ratio * e.z);
		r_postproc = min((float)(ratio * e.r), (float)e.r * .5f);
	}
#define F(func) func(z_postproc, r_postproc, color);
	F(DrawCircleLinesV);
	F(DrawCircleV);
#undef F
}

/// <summary>
/// Universal gravitational constant (units: LLL/T/T/M).
/// </summary>
constexpr double G = .1;

/// <summary>
/// Time step (T per frame).
/// </summary>
constexpr double DT = 0.005;

/// <summary>
/// Force on the left particle (l) due to the right particle (r).
/// </summary>
/// <returns>Force (units: ML/T/T/T)</returns>
static C newton_gravity(Dyn::Entry const& l, Dyn::Entry const& r)
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
		constexpr int M = 25;
		// Force per mass (integrated).
		C fpm;
		// Number of pieces sampled in the left crescent (one-sided lune---just "lune").
		int n{};
		// Circular intersection.
		CircularIntersection sect(l.z, l.r, r.z, r.r);
		// Computation of lunar force.
		auto infinitesimal = [&](C const& p)
			{
				// `p` is sampled from a certain square in
				// a reoriented coordinate system, where
				// the left circle is centered at the origin, and
				// the right circle as at (`as`, 0). Lengths have
				// not changed.

				if (!sect.left(p) || sect.right(p)) return;
				n++;
				C arm = as - p; double dist = abs(arm);
				// [***] `fpm` will be missing the factors of: M (= # trials), G, dm.
				// [***] `dm` is as yet unavailable; it's computed soon.
				fpm += 1 / dist / dist / dist / M * sect.unrotate(arm);
			};
		// Boom.
		for (int i = M - 1; i >= 0; i--)
			sect.monte(hh.next(), infinitesimal);
		// If no sample hit, the integration failed.
		if (!n) return 0;
		// [***] Compute dm by the ratio -> dm : m = 1 : n
		// (where dm: infinitesimal mass, m: mass of left particle,
		// n: number of samples hit).
		// You can see why I have left out the computation of dm:
		// the number of samples hit, n, is unavailable while "n++".
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

/// <summary>
/// Compute the largest absolute value between the
/// respective differences of the real and imaginary
/// components of the given complex numbers `a` and `b`.
/// 
/// Complex Largest Absolute Deviation.
/// </summary>
static double clad(C const& a, C const& b)
{
	return std::max(
		abs(a.real() - b.real()),
		abs(a.imag() - b.imag())
	);
}

/// <summary>
/// Judge the two calculated position values that should ideally be
/// identical (but would be different if the system was too violent).
/// 
/// The specification is in the `Dyn::Driver` structure documentation.
/// Basically, +1 expresses judgement of safety (so use a larger time
/// step); -1 expresses concern (use a finer time step and retry the computation
/// as appropriate); 0 expresses neutrality.
/// 
/// Strong vs. weak: this distinction is explained in the
/// `beasons` namespace docs. See Beasons.h.
/// Beason's method, by the way, is the chosen method of integration.
/// 
/// Basically, the strong one is the one that will be
/// substitued in for the position value of the particle at the next time step,
/// and the weak one is a duplicate calculation that is
/// only provided for the estimation of error.
/// Again, in the absense of error, strong should nearly equal weak.
/// </summary>
static int judge_z(C const& strong, C const& weak)
{
	double c = clad(strong, weak);
	// Units: L.
	if (c > 0.001) return -1; // try finer time step.
	else if (c < 0.0001) return +1; // suggest coarser time step.
	else return 0;
}

/// <summary>
/// Like `judge_z`, judge the velocities.
/// </summary>
static int judge_v(C const& strong, C const& weak)
{
	double c = clad(strong, weak);
	// Units: L/T.
	if (c > 0.001) return -1; // try finer time step.
	else if (c < 0.0001) return +1; // suggest coarser time step.
	else return 0;
}

static Dyn make()
{
	Dyn dyn;
	dyn.par.dt = DT;
	dyn.drv.judge_z = judge_z;
	dyn.drv.judge_v = judge_v;
	{
		// Generate this many (n) particles.
		int constexpr n = 125;
		auto seed = []() { std::random_device dev; return dev(); }();
		auto rng = std::mt19937(seed);
		C const rot = std::polar(1., PI64 / 3);
		for (int i = n - 1; i >= 0; i--)
		{
			std::uniform_real_distribution<> v(-10, 10), r(1.0, 5.0);
			std::cauchy_distribution<> z(0., 30.), m(20, 7.); // center; scale.
			auto sq = [](double a) { return a * a; };
#define sca(d) d(rng)
#define vec(d) C(sca(d), sca(d))
			Dyn::Entry e;
			e.z = vec(z), e.v = vec(v) + rot / abs(e.z) * e.z, e.a = 0;
			e.m = sq(sca(m)) + 1., e.r = sq(sca(r)) + 1.;
#undef vec
#undef sca
			dyn.tab.push_back(e);
		}
		dyn.drv.pair_force = newton_gravity;
		// It is here where all accelerations are computed
		// for before the first iteration, and where the
		// the total mass (dyn.m_mass) is computed.
		dyn.precompute();
	}
	return dyn;
}

static Dyn make_set1()
{
	Dyn dyn;
	dyn.par.dt = DT;
	dyn.drv.judge_z = judge_z;
	dyn.drv.judge_v = judge_v;
	Dyn::Entry e0, e1;
	e0.z = -10., e1.z = -e0.z;
	e0.m = 30., e1.m = e0.m;
	e0.r = 10., e1.r = e0.r;
	Dyn::Entry e2(e1);
	e2.z = 20.;
	e2.r /= 4;
	dyn.tab.push_back(e0);
	dyn.tab.push_back(e1);
	dyn.tab.push_back(e2);
	dyn.drv.pair_force = newton_gravity;
	dyn.precompute();
	return dyn;
}

static double kinetic_energy(Dyn const& dyn)
{
	double ke{};
	for (int i = dyn.n() - 1; i >= 0; i--)
		ke += std::norm(dyn[i].v) * dyn[i].m;
	return ke / 2;
}

static void universal_force(Dyn& dyn)
{
	for (int i = dyn.n() - 1; i >= 0; i--)
	{
		auto& e = dyn[i];
		// Unphysical effect(s).

		// "Drag"
		//double av = abs(e.v);
		//double av2 = 350. * tanh(av / 350.);
		//e.v *= av2 / av;
	}
}

int wWinMain(void* _0, void* _1, void* _2, int _3)
{
	// auto sim = make_set1;
	auto sim = make;

	// Simulation (dyn)
	Dyn dyn = sim();

	// Rendering
	float constexpr px_per_l = 1.f;

	// Misc.
	int constexpr reset_at_sec = 180;
	int resets = 0;
	double last_reset_s = 0;

	// Raylib.
	constexpr int fps_target = 60;
	InitWindow(600, 600, "Gravity");
	SetTargetFPS(fps_target);

	// If there's time, schedule more calls to the simulation per frame.
	auto const load_good = [=]() { return (double)GetFPS() >= 0.90 * fps_target; };
	auto const load_terrible = [=]() { return (double)GetFPS() <= 0.65 * fps_target; };
	int const scheduling_levelup_at = 20;
	int const scheduling_level_cap_excl = 15;
	int scheduling_mood = 0;
	auto const calls_per_frame = [&]() { return 1 + scheduling_mood / scheduling_levelup_at; };
	auto const up_mood = [&]()
		{
			scheduling_mood++;
			if (calls_per_frame() >= scheduling_level_cap_excl)
				scheduling_mood--;
		};
	auto const down_mood = [&]()
		{
			if (--scheduling_mood < 0) scheduling_mood++;
		};

	while (!WindowShouldClose())
	{
		if (IsKeyPressed(KEY_R))
		{
			// reset simulation
			dyn = sim();
			last_reset_s = GetTime();
			resets = 0;
			scheduling_mood = 0;
		}
		else
		{
			// Regularly reset.
			double time = GetTime() - last_reset_s;
			int quo = (int)(time / reset_at_sec);
			if (quo > resets)
			{
				dyn = sim();
				scheduling_mood = 0;
			}
			resets = std::max(quo, resets);
		}

		for (int i = calls_per_frame() - 1; i >= 0; i--)
		{
			dyn.step();
			dyn.bias();
			universal_force(dyn);
		}

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
			auto ke = kinetic_energy(dyn);
			char msg_ke[500];
			snprintf(msg_ke, sizeof(msg_ke),
				"KE: %.4G MLL/T/T\n"
				"dt: %.6f T/try\n"
				"tries per frame: %d",
				ke, dyn.par.dt, calls_per_frame()
			);
			DrawText(msg_ke, 16, 40, 20, BLACK); // x, y, font size (px)
		}
		EndDrawing();

		// Reflect upon the performance.
		if (load_good()) up_mood();
		else if (load_terrible()) down_mood();
	}

	CloseWindow();
	return 0;
}
