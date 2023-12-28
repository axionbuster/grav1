#pragma once

#include "Include.h"

inline Cf downgrade(C c)
{
	return Cf((float)c.real(), (float)c.imag());
}

inline Vector2 c2vec(Cf c)
{
	return Vector2{ c.real(),c.imag() };
}

inline Vector2 c2vec(C c)
{
	return c2vec(downgrade(c));
}

inline Cf vec2c(Vector2 v)
{
	return Cf(v.x, v.y);
}

inline Cf round(Cf c)
{
	return Cf(round(c.real()), round(c.imag()));
}

/// <summary>
/// Some parameters for visualization.
/// </summary>
struct VisParam
{
	/// <summary>
	/// Scale factor (pixels : L unit), where
	/// an "L unit" is the internal (physics engine) length unit.
	/// </summary>
	float sc_px_per_l{ 1.0f };
	float radius_px{ 5.0f };
	float arrow_max_px{ 30.0f };
	float arrow_tip_len_ratio{ 1.0f / 3 };
	Color circle_color{ BLACK };
	Color label_color{ BLACK };
};

/// <summary>
/// Interactive visualizer; primitives and shared data to put stuff onto the screen.
/// </summary>
class Vis
{
public:
	VisParam param{};
	/// <summary>
	/// In the screen-coordinate system, where the origin of the internal
	/// coordinate system should be located, measured in pixels.
	/// 
	/// Safe to use the midpoint of the screen width and height dimensions.
	/// 
	/// Intended for frequent change.
	/// </summary>
	Cf origin_px{};

	Vis() {}
	Vis(VisParam param) : param(param) {}
	Vis(VisParam param, Cf origin_px) : param(param), origin_px(origin_px) {}

	/// <summary>
	/// Scale and translate a point in the internal coordinate system
	/// in L units to the screen coordinate system in pixels.
	/// </summary>
	/// <param name=""></param>
	Cf locate(Cf z) const
	{
		return z * param.sc_px_per_l + origin_px;
	}

	/// <summary>
	/// Plot a circle centered at z in the internal coordinate system
	/// in L units with a radius that is the product of rad_mult and
	/// param.radius_px.
	/// </summary>
	/// <param name=""></param>
	/// <param name="rad_mult"></param>
	void plot(Cf z, float rad_mult = 1.0f) const
	{
		DrawCircleLinesV(c2vec(locate(z)), param.radius_px * rad_mult, param.circle_color);
	}

	void plot(C z, double rad_mult = 1.0f) const { plot(downgrade(z), (float)rad_mult); }

	void label(Cf z, const char* str) const
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

	void label(C z, const char* str) const { label(downgrade(z), str); }

	/// <summary>
	/// Plot an "arrow" beginning at z with the direction and magnitude of v
	/// in "internal" (L) units.
	/// </summary>
	/// <param name=""></param>
	/// <param name=""></param>
	void arrow_at(Cf z, Cf v) const
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

	void arrow_at(C z, C v) const { arrow_at(downgrade(z), downgrade(v)); }
};
