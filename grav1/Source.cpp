#include "Dyn.h"
#include "Vis.h"

#include <functional>
#include <random>

constexpr auto FPS = 60;

/// <summary>
/// Draw things to the screen and talk to the human,
/// while advancing the simulation state.
/// </summary>
/// <param name="">Physical table (particles and states);
/// accelerations must have been computed.</param>
/// <param name="">Simulation parameters, which this function might change.</param>
/// <param name="">Visualizer object, which this function might touch.</param>
/// <returns>Suggested return value for the process (e.g., 0 means success).</returns>
int window(Dyn&, Param&, Vis&);

/// <summary>
/// Generate a random number (stored in the heap).
/// 
/// This structure is Callable.
/// </summary>
template<typename State = std::mt19937>
struct HeapRandom
{
	/// <summary>
	/// Given a seed (specified by the internal state structure),
	/// initialize a random state.
	/// 
	/// Safe to call a std::random_device instance for the seed.
	/// </summary>
	/// <typeparam name="Seed"></typeparam>
	/// <param name="s"></param>
	template<typename Seed>
	HeapRandom(Seed s) : mt(std::make_shared<State>(s)) {}
	/// <summary>
	/// Generate a random number (double)
	/// from a given distribution (dist : Dist).
	/// </summary>
	/// <typeparam name="Dist">A type for the distribution.</typeparam>
	/// <param name="dist">The distribution.</param>
	/// <returns></returns>
	template<typename Dist>
	double operator()(Dist& dist) { return dist(*mt); }
private:
	/// <summary>
	/// Internal random number generator state.
	/// 
	/// (Mersenne Twister).
	/// </summary>
	std::shared_ptr<State> mt;
};

int wWinMain(void* _h, void* _h2, void* _c, int _s)
{
	int n = 700;
	Param param{};
	Vis vis;
	auto tab(new DynEntry[n]{});

	Dyn d(n, param, std::shared_ptr<DynEntry[]>(tab));

	// Some policies
	auto test_replace = [](DynEntry const& d)
		{
#define cfinite(c) (std::isfinite((c).real()) && std::isfinite((c).imag()))
			return
				!(cfinite(d.z) && cfinite(d.v) && cfinite(d.a)
					&& std::isnormal(d.m) && abs(d.z) < 10000);
#undef cfinite
		};
	std::random_device rd;
	auto replace0 = [](HeapRandom<> hr)
		{
			typedef std::uniform_real_distribution<> urdist;
			urdist posdist(-300, 300), massdist(0.4, 0.9), veldist(0, 1),
				vel_theta_dist(0.0, PI / 3);
			DynEntry d;
			d.z = C(hr(posdist), hr(posdist));
			d.m = hr(massdist);
			d.v = std::polar(hr(veldist) / abs(d.z), hr(vel_theta_dist)) * d.z;
			d.ljsigma = std::sqrt(d.m) * 10.0;
			return d;
		};
	auto replace = std::bind(replace0, HeapRandom<>(rd()));
	d.driver.test_replace = test_replace;
	d.driver.replace = replace;

	// Integration method requires calculation of accelerations
	// before the first iteration.
	d.accelall();

	return window(d, param, vis);
}

int window(Dyn& dyn, Param& sim_param, Vis& vis)
{
	// Raylib.
	InitWindow(800, 600, "a");
	SetTargetFPS(FPS);

	uint32_t fr = 0;
	while (!WindowShouldClose())
	{
		// Do the calculations.
		dyn.iterall();
		dyn.center_vz();

		int n = dyn.n;
		auto const& tab = dyn.tab;

		// Draw.

		vis.origin_px = Cf(
			GetScreenWidth() / 2.0f, GetScreenHeight() / 2.0f
		);

		BeginDrawing();
		{
			ClearBackground(WHITE);

			for (int i = n - 1; i >= 0; i--)
			{
				vis.plot(tab[i].z, tab[i].m);
				vis.arrow_at(tab[i].z, tab[i].v);
			}

			// where (px): x, y.
			DrawFPS(16, 16);
		}
		EndDrawing();

		fr++;
	}
	CloseWindow();
	return 0;
}
