#include "Header.h"
#include "Vis.h"

#include <random>

constexpr auto FPS = 60;

int window(Dyn&, Param const&, Vis&);

int wWinMain(void* _h, void* _h2, void* _c, int _s)
{
	int n = 500;
	Param param{};
	Vis vis;
	auto kin(new KinEntry[n]);

	{
		// Lots of random tweakings here and there...
		std::random_device rd;
		std::mt19937 mt(rd());
		std::normal_distribution<double> dist;
#define d() dist(mt)
#define cd() C(d(), d())
		for (int i = n - 1; i >= 0; i--)
			kin[i].m = abs(d()) + 0.01, kin[i].z = cd() * 25.0;
		// Cool "vortex" effect.
		C curl0(std::polar(abs(d()) * 40.0, d()));
		for (int i = n - 1; i >= 0; i--)
			kin[i].v = curl0 * kin[i].z / abs(kin[i].z);
#undef cd
#undef d
	}
	// Strictly improper, but my debugger can't inspect the contents
	// of kin here without me exposing the inner object (kin as a raw array).
	Dyn d(n, param, std::shared_ptr<KinEntry[]>(kin));

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
		dyn.iterall();

		if (fr % (FPS * 5) == 0)
			dyn.reset_centroid();

		int n = dyn.n;
		auto const& kin = dyn.kin;

		vis.origin_px = Cf(
			GetScreenWidth() / 2.0f, GetScreenHeight() / 2.0f
		);

		BeginDrawing();
		{
			ClearBackground(WHITE);

			for (int i = n - 1; i >= 0; i--)
			{
				vis.plot(kin[i].z, kin[i].m);
				vis.arrow_at(kin[i].z, kin[i].v);

				// Numerical label above particle.
				char label[33]{};
				if (_itoa_s(i, label, 10) != 0)
					label[0] = '?', label[1] = '\0';
				vis.label(kin[i].z, label);
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
	// Newton's gravity has an inherent singularity at r = 0 (zero distance).
	// I have to do something sensible in that situation.
	static auto guard0 = [=](double x) { return x < param.guard0_dist ? param.guard0_dist : x; };
	static auto invsq = [=](C z) {
		double w = abs(z);
		// Yes, divide three times.
		// First to normalize the vector z, creating z / |z|
		// (but we only return the scalar multiplier part, note!)
		// then divide twice more so we have (z / |z|) / |z| / |z|,
		// which is by algebra unit(z) / |z|^2.
		// ... Inverse square law, as desired.
		return guard0(1 / w / w / w); // type: (double)
		};
	C a;
	for (int j = n - 1; j >= 0; j--)
	{
		if (i == j) continue;
		// Newton's law of gravitation.
		C zj = kin[j].z, r = zj - z, aij = param.g * kin[j].m * invsq(r) * r;
		a += aij;
	}
	return a;
}

void Dyn::accelall() const
{
	for (int i = n - 1; i >= 0; i--)
		kin[i].a = accel(i, kin[i].z);
}

KinEntry Dyn::iter(int i) const
{
	static auto speed_limit = [=](C v) {
		double w = abs(v);
		if (w == 0.0) return v;
		// Take this ratio, which should be close to 0 if v is small.
		double ww = w / param.speed_limit;
		// If small enough, don't do anything.
		// If too large, cap it to a certain proportion of the speed limit.
		return ww < 0.99 ? v : 0.99 / ww * v;
		};

	// Method of leapfrog integration: conserves energy.
	double dt = param.dt, dthalf = dt / 2;
	KinEntry k(kin[i]);
	C v1_(k.v + k.a * dthalf), v1(speed_limit(v1_)), z2(k.z + v1 * dt);
	C a2(accel(i, z2));
	C v2_(v1 + a2 * dthalf), v2(speed_limit(v2_));
	k.z = z2, k.v = v2, k.a = a2;
	return k;
}

void Dyn::iterall() const
{
	for (int i = n - 1; i >= 0; i--)
		kin[i] = iter(i);
}

C Dyn::centroid() const
{
	// Welford's online algorithm.
	// Difficult to parallelize, but at least numerically stable.
	C welford;
	for (int i = n - 1, m = 0; i >= 0; i--)
		welford += (kin[i].z - welford) / (double)++m;
	return welford;
}

void Dyn::reset_centroid() const
{
	C c(centroid());
	for (int i = n - 1; i >= 0; i--) kin[i].z -= c;
}
