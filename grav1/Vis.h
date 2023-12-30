#pragma once

#include "Include.h"

/// <summary>
/// Downgrade the precision of the floating point from 64 bit to 32 bit.
/// </summary>
/// <param name="c"></param>
Cf downgrade(C c);

/// <summary>
/// Some parameters for visualization.
/// </summary>
struct VisParam
{
	/// <summary>
	/// Scale factor (pixels : L unit), where
	/// an "L unit" is the internal (physics engine) length unit.
	/// </summary>
	float sc_px_per_l{ 0.50f };
	float radius_px{ 15.0f };
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
	Cf locate(Cf z) const;

	/// <summary>
	/// Plot a circle centered at z in the internal coordinate system
	/// in L units with a radius that is the product of rad_mult and
	/// param.radius_px.
	/// </summary>
	/// <param name=""></param>
	/// <param name="rad_mult"></param>
	void plot(Cf z, float rad_mult = 1.0f) const;

	void plot(C z, double rad_mult = 1.0f) const { plot(downgrade(z), (float)rad_mult); }

	/// <summary>
	/// Put a label over a position in internal (L) units / coordinates.
	/// </summary>
	/// <param name="z">where (L units / coordinates)</param>
	/// <param name="str"></param>
	void label(Cf z, const char* str) const;

	void label(C z, const char* str) const { label(downgrade(z), str); }

	/// <summary>
	/// Plot an "arrow" beginning at z with the direction and magnitude of v
	/// in "internal" (L) units.
	/// </summary>
	/// <param name=""></param>
	/// <param name=""></param>
	void arrow_at(Cf z, Cf v) const;

	void arrow_at(C z, C v) const { arrow_at(downgrade(z), downgrade(v)); }
};
