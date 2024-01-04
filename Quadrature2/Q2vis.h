#pragma once
#include "Header.h"

/// <summary>
/// Visualization logic
/// </summary>
namespace vis {
	/// <summary>
	/// Drop the precision.
	/// </summary>
	/// <param name="c"></param>
	/// <returns></returns>
	inline Cf downgrade(C c) { return Cf((float)c.real(), (float)c.imag()); }
	/// <summary>
	/// Convert a 32-bit complex number to `Vector2` (float) from Raylib.
	/// </summary>
	/// <param name="c"></param>
	/// <returns></returns>
	inline Vector2 v2(Cf c) { return Vector2{ c.real(), c.imag() }; }
	/// <summary>
	/// Convert a 64-bit complex number to `Vector2` (float) from Raylib.
	/// </summary>
	/// <param name="c"></param>
	/// <returns></returns>
	inline Vector2 v2(C c) { return v2(downgrade(c)); }
	/// <summary>
	/// Plot a single point in screen coordinates (pixels).
	/// </summary>
	/// <param name="v">Point in screen coordinates, measured in pixels</param>
	/// <param name="c">Color</param>
	inline void plot(Cf v, Color c = BLACK)
	{
		Cf const dim(3.f, 3.f); // pixels (side length)
		v -= dim * 0.5f;
		Vector2 v2{ v.real(), v.imag() },
			dim2{ dim.real(), dim.imag() };
		DrawRectangleV(v2, dim2, c);
	}
};
