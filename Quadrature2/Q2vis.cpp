#include "Q2vis.h"

void vis::plot(Cf v, Color c)
{
	Cf const dim(3.f, 3.f); // pixels (side length)
	v -= dim * 0.5f;
	Vector2 v2{ v.real(), v.imag() },
		dim2{ dim.real(), dim.imag() };
	DrawRectangleV(v2, dim2, c);
}

void vis::circle(Cf o, float r, Color const* fill, Color const* border)
{
	Vector2 v{ o.real(), o.imag() };
	if (fill)
		DrawCircleV(v, r, *fill);
	if (border)
		DrawCircleLinesV(v, r, *border);
}
