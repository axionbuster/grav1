#include "Vis.h"

Cf downgrade(C c)
{
	return Cf((float)c.real(), (float)c.imag());
}

inline static Vector2 c2vec(Cf c)
{
	return Vector2{ c.real(), c.imag() };
}

inline static Vector2 c2vec(C c)
{
	return c2vec(downgrade(c));
}

inline static Cf vec2c(Vector2 v)
{
	return Cf(v.x, v.y);
}

Cf Vis::locate(Cf z) const
{
	return z * param.sc_px_per_l + origin_px;
}

void Vis::plot(Cf z, float rad_mult) const
{
	DrawCircleLinesV(c2vec(locate(z)),
		param.radius_px * rad_mult,
		param.circle_color);
}

inline static Cf round(Cf c)
{
	return Cf(round(c.real()), round(c.imag()));
}

void Vis::label(Cf z, const char* str) const
{
	// Some boilerplate.
	Font f = GetFontDefault();
	float font_size = 12.0f;
	float spacing = 1.0f;
	Color tint = param.label_color;

	// Center and push upward.
	Vector2 dim = MeasureTextEx(f, str, font_size, spacing);
	Cf at_px = locate(z) - Cf(dim.x, dim.y) * 0.5f - 16.0if;
	Cf at_px2 = round(at_px);

	DrawTextEx(f, str, c2vec(at_px2), font_size, spacing, tint);
}

void Vis::arrow_at(Cf z, Cf v) const
{
	static const Cf ang = std::polar(1.0f, 30.0f * PI / 180),
		ang2 = std::conj(ang);
	// In pixel units, cap the length of the arrow.
	Cf vpx = v * param.sc_px_per_l;
	float apx = abs(vpx);
	if (apx >= param.arrow_max_px)
	{
		vpx *= param.arrow_max_px / apx;
		apx = param.arrow_max_px;
	}
	// In pixel units, construct both segments
	// that make up the "tip" part of the arrow.
	Cf vpx_scaled_down = param.arrow_tip_len_ratio * vpx;
	Cf tip1 = vpx_scaled_down * ang, tip2 = vpx_scaled_down * ang2;
	// Compute their locations on the screen.
	Cf z1 = locate(z), v1 = z1 + vpx,
		tip3 = v1 - tip1, tip4 = v1 - tip2;
	DrawLineEx(c2vec(z1), c2vec(v1), 1.0f, BLACK);
	DrawLineEx(c2vec(v1), c2vec(tip3), 1.0f, BLACK);
	DrawLineEx(c2vec(v1), c2vec(tip4), 1.0f, BLACK);
}
