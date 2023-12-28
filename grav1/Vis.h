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

struct VisParam
{
	float sc_px_per_l{};
};

inline void plot(VisParam param, Cf z_l, Cf o_px) 
{
	Cf z_o_px = z_l * param.sc_px_per_l + o_px;
	// radius: 5.0 pixels.
	DrawCircleLinesV(c2vec(z_o_px), 5.0f, BLACK);
}
