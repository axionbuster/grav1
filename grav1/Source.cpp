#include "Header.h"
#include "Vis.h"

#include <random>

constexpr auto FPS = 60;

int window(Dyn&, Param const&, Vis&);

int wWinMain(void* _h, void* _h2, void* _c, int _s)
{
	int n = 1000;
	Param param{};
	Vis vis;
	auto tab(new DynEntry[n]{});

	{
		// Lots of random tweakings here and there...
		std::random_device rd;
		std::mt19937 mt(rd());

		// units of measurement in this block:
		// "internal" units (not display units).

		// orientation: left-handed.
		const double width = 500.0,
			height = width,
			left = -width / 2, right = width / 2,
			top = -height / 2, bottom = height / 2;
		const int rows = 20,
			columns = 20;

#define d(dist) dist(mt)
#define urdist std::uniform_real_distribution<double>
#define make(what, dist0, dist1) what.z = C(dist0(mt), dist1(mt))
		int cells = rows * columns,
			quo = n / cells,
			i = n - 1;
		urdist mass(1.0, 5.0);
		for (int r = rows - 1; r >= 0; r--)
		{
			const double row_top = top + height / rows * r,
				row_bottom = top + height / rows * (r + 1);
			urdist distr(row_top, row_bottom);
			for (int c = columns - 1; c >= 0; c--)
			{
				const double column_left = left + width / columns * c,
					column_right = left + width / columns * (c + 1);
				urdist distc(column_left, column_right);
				for (int j = quo - 1; j >= 0; j--)
				{
					make(tab[i], distc, distr);
					tab[i].m = mass(mt);
					tab[i].ljsigma = tab[i].m * 5.0;
					i--;
				}
			}
		}
		urdist distr(top, bottom), distc(left, right);
		while (i >= 0)
		{
			make(tab[i], distc, distr);
			tab[i].m = mass(mt);
			tab[i].ljsigma = tab[i].m * 5.0;
			i--;
		}

		// Some masses ... they are made extra massive and reactive.
		const double prob = 0.01;
		urdist bernoulli;
		for (i = n - 1; i >= 0; i--)
			if (bernoulli(mt) < prob)
				tab[i].m *= 10.0, tab[i].ljsigma *= 10.0;
#undef make
#undef urdist
#undef d
	}
	Dyn d(n, param, std::shared_ptr<DynEntry[]>(tab));

	// Integration method requires calculation of accelerations
	// before the first iteration.
	d.accelall();

	return window(d, param, vis);
}

int window(Dyn& dyn, Param const& sim_param, Vis& vis)
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

C Dyn::accel(int i, C z) const
{
	C accel;
	for (int j = n - 1; j >= 0; j--)
	{
		if (i == j) continue;
		C zj = tab[j].z, r = zj - z;
		// Newton's gravity has an inherent singularity at r = 0 (zero distance).
		// I have to do something sensible in that situation.
		double w = std::max(abs(r), param.guard0_dist);
		// Newton's law of gravity.
		C aij_newtgrav = param.g * tab[j].m * (1 / w / w / w) * r;
		// Lennard-Jones potential (gradient).
		// Let h = signa(j) / |r|.
		// where j is the "other" particle.
		// Note that h is dimensionless, since both |r| and sigma are lengths.
		// Now, compute the 7th and 13th powers of h.
		double h = tab[j].ljsigma / w,
			h3 = h * h * h, h6 = h3 * h3, h12 = h6 * h6;
		C aij_lj = 4 * param.lj_force_unit * (h * h6 - h * h12) / tab[i].m * r;
		accel += aij_newtgrav + aij_lj;
	}
	return accel;
}

void Dyn::accelall() const
{
	std::unique_ptr<DynEntry[]> copy(new DynEntry[n]{});
	for (int i = n - 1; i >= 0; i--)
	{
		copy[i] = tab[i];
		copy[i].a = accel(i, copy[i].z);
	}
	for (int i = n - 1; i >= 0; i--)
		tab[i] = copy[i];
}

DynEntry Dyn::iter(int i) const
{
	// Limit the magnitude of a vector.
	auto limit = [=](double lim, C x) {
		double w = abs(x);
		if (w == 0.0) return x;
		// Take this ratio, which should be close to 0 if x is small.
		double ww = w / lim,
			mag = lim / w * (2.0 / PI) * atan(ww);
		return mag * x;
		};
	auto speed_limit = [=](C v) { return limit(param.speed_limit, v); };
	auto accel_limit = [=](C a) { return limit(param.accel_limit, a); };

	// Method of leapfrog integration: conserves energy.
	double dt = param.dt;
	DynEntry t(tab[i]);
	C v1_(t.v + dt / 2 * t.a), v1(speed_limit(v1_)), z2(t.z + v1 * dt);
	C a2_(accel(i, z2)), a2(accel_limit(a2_));
	C v2_(v1 + dt / 2 * a2), v2(speed_limit(v2_));
	t.z = z2, t.v = v2, t.a = a2;
	return t;
}

void Dyn::iterall() const
{
	std::unique_ptr<DynEntry[]> copy(new DynEntry[n]{});
	for (int i = n - 1; i >= 0; i--)
		copy[i] = iter(i);
	for (int i = n - 1; i >= 0; i--)
		tab[i] = copy[i];
}

void Dyn::center_vz() const
{
	// Welford's online algorithm.
	// Difficult to parallelize, but at least numerically stable.
	C barycenter, velbias;
	for (int i = n - 1, m = 0; i >= 0; i--)
	{
		double o = ++m;
		barycenter += (tab[i].z * tab[i].m - barycenter) / o;
		velbias += (tab[i].v - velbias) / o;
	}
	barycenter /= mass;
	for (int i = n - 1; i >= 0; i--) tab[i].z -= barycenter, tab[i].v -= velbias;
}
