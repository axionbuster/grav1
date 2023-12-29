#include "Header.h"
#include "Vis.h"

#include <random>

constexpr auto FPS = 60;

int window(Dyn&, Param const&, Vis&);

int wWinMain(void* _h, void* _h2, void* _c, int _s)
{
	int n_half = 200, n = n_half * 2;
	Param param{};
	Vis vis;
	auto kin(new KinEntry[n]{});

	{
		// Lots of random tweakings here and there...
		std::random_device rd;
		std::mt19937 mt(rd());
		std::normal_distribution<double> dist;
#define d() dist(mt)
#define cd() C(d(), d())
		double largest_mass{};
		for (int i = n - 1; i >= 2; i--)
		{
			double m = abs(d()) * 0.2 + 0.1;
			largest_mass = std::max(largest_mass, m);
			kin[i].m = m, kin[i].z = cd() * 25.0,
				kin[i].v = cd() * 0.1;
		}
		kin[0].m = largest_mass + 10.0;
		kin[1].m = largest_mass * 0.1 + 10.0;
		C second_center(cd() * 300.0);
		C second_centers_velocity(cd() * 5.0);
		kin[1].z += second_center;
		kin[1].v += second_centers_velocity;
		for (int i = n - 1; i >= n_half; i--)
			kin[i].z += second_center;
#undef cd
#undef d
	}
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
		// Do the calculations.
		dyn.iterall();
		dyn.reset_bary();

		int n = dyn.n;
		auto const& kin = dyn.kin;

		// Draw.

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
				char label[33];
				if (_itoa_s(i, label, 10) != 0)
					label[0] = '?', label[1] = '\0';
				// vis.label(kin[i].z, label);
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
	auto cbrecip_guard0 = [=](double x) {
		x = std::max(x, param.guard0_dist);
		return 1 / x / x / x;
		};
	C specific_force;
	for (int j = n - 1; j >= 0; j--)
	{
		if (i == j) continue;
		// Newton's law of gravitation.
		C zj = kin[j].z, r = zj - z;
		double w = cbrecip_guard0(abs(r));
		C aij = param.g * kin[j].m * w * r;
		specific_force += aij;
	}
	return specific_force;
}

void Dyn::accelall() const
{
	std::unique_ptr<KinEntry[]> copy(new KinEntry[n]);
	for (int i = n - 1; i >= 0; i--)
		copy[i].a = accel(i, (copy[i] = kin[i], copy[i].z));
	for (int i = n - 1; i >= 0; i--)
		kin[i] = copy[i];
}

KinEntry Dyn::iter(int i) const
{
	auto speed_limit = [=](C v) {
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
	std::unique_ptr<KinEntry[]> copy(new KinEntry[n]);
	for (int i = n - 1; i >= 0; i--)
		copy[i] = iter(i);
	for (int i = n - 1; i >= 0; i--)
		kin[i] = copy[i];
}

C Dyn::bary() const
{
	// Welford's online algorithm.
	// Difficult to parallelize, but at least numerically stable.
	C welford;
	for (int i = n - 1, m = 0; i >= 0; i--)
		welford += (kin[i].z * kin[i].m - welford) / (double)++m;
	return welford;
}

void Dyn::reset_bary() const
{
	C c(bary());
	for (int i = n - 1; i >= 0; i--) kin[i].z -= c;
}
