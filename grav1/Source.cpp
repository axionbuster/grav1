#include "Header.h"
#include "Vis.h"

int window(Dyn&, Param const&, VisParam const&);

int wWinMain(void* _h, void* _h2, void* _c, int _s)
{
	int n = 2;
	Param p0{};
	VisParam p1{};
	p1.sc_px_per_l = 10.0f;
	auto kin(new KinEntry[2]);
	kin[0].m = 1.0, kin[0].v = 0.0 + 0.0i, kin[0].z = -3.1 - 0.3i;
	kin[1].m = 0.5, kin[1].v = 0.0 + 0.0i, kin[1].z = 0.0 + 0.0i;
	// Strictly improper, but my debugger can't inspect the contents
	// of kin here without me exposing the inner object.
	Dyn d(2, p0, std::shared_ptr<KinEntry[]>(kin));

	d.accelall();

	return window(d, p0, p1);
}

int window(Dyn& dyn, Param const& sim_param, VisParam const& vis_param)
{
	InitWindow(800, 600, "a");
	SetTargetFPS(60);
	while (!WindowShouldClose())
	{
		dyn.iterall();

		int n = dyn.n;
		auto const& kin = dyn.kin;

		Cf screen_midpoint_px(
			GetScreenWidth() / 2.0f, GetScreenWidth() / 2.0f
		);

		BeginDrawing();
		{
			ClearBackground(WHITE);

			for (int i = n - 1; i >= 0; i--)
			{
				Cf z = downgrade(kin[i].z);
				plot(vis_param, z, screen_midpoint_px);
			}
		}
		EndDrawing();
	}
	CloseWindow();
	return 0;
}

C Dyn::accel(int i, C z) const
{
	static auto guard0 = [=](double x) { return x < param.guard0_dist ? param.guard0_dist : x; };
	static auto invsq = [=](C z) {
		double w = abs(z);
		// Yes, divide three times.
		// First two to normalize the vector z, creating z / |z|
		// (but we only return the scalar multiplier part, note!)
		// then divide once more so we have (z / |z|) / |z|,
		// which is by algebra unit(z) / |z|^2.
		// ... Inverse square law, as desired.
		return guard0(1 / w / w / w); // double
		};
	C a;
	for (int j = n - 1; j >= 0; j--)
	{
		if (i == j) continue;
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
	double dt = param.dt, dthalf = dt / 2;
	KinEntry k(kin[i]);
	C v1(k.v + k.a * dthalf), z2(k.z + v1 * dt);
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
