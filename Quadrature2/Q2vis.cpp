#include "Q2vis.h"

void vis::plot(Cf v, Color c)
{
	Cf const dim(5.f, 5.f);
	v -= dim * 0.5f;
	Vector2 v2{ v.real(), v.imag() },
		dim2{ dim.real(), dim.imag() };
	DrawRectangleV(v2, dim2, c);
}

static void line(Cf from, Cf to, float thick = 1.f, Color c = BLACK)
{
	Vector2 from2{ from.real(), from.imag() },
		to2{ to.real(), to.imag() };
	DrawLineEx(from2, to2, thick, c);
}

static void axis1(Cf v, Cf interval, int n, char const* format, double show_unit = 0.f)
{
	Cf const tick = 5.if;
	Cf const rot = interval / abs(interval);
	Cf const neg_extr = -interval * (float)n,
		pos_extr = interval * (float)n;
	line(v + neg_extr, v + pos_extr);

	for (float i = (float)-n; i <= n; i++)
	{
		if (i == 0) continue;
		Cf const at = interval * i,
			to1 = at + rot * tick * .5f,
			to2 = at - rot * tick * .5f;
		line(v + to1, v + to2);

		if (show_unit == 0) continue;
		// Some font stuff...
		Font font = GetFontDefault();
		float font_size =
			(i == 1 || i == -1) ? 16.f : 10.f;
		float spacing = 1.0f;
		Color tint = BLACK;

		// Plan printing.
		char text[256];
		snprintf(text, sizeof(text), format, i * show_unit);
		Vector2 text_dim = MeasureTextEx(font, text, font_size, spacing);
		Cf text_dim_cf(text_dim.x, text_dim.y);

		// Center.
		Cf text_loc = at - text_dim_cf * .5f,
			text_loc_v = v + text_loc;
		Vector2 text_loc_v2{ text_loc_v.real(), text_loc_v.imag() };

		// Print.
		DrawTextEx(font, text, text_loc_v2, font_size, spacing, tint);
	}
}

void vis::axes(Cf v, Cf re, Cf im, int n, double label)
{
	axis1(v, re, n, "%.2f", label);
	axis1(v, im, n, "%.2fi", label);
}

void vis::circle(Cf o, float r, Color const* fill, Color const* border)
{
	Vector2 v{ o.real(), o.imag() };
	if (fill)
		DrawCircleV(v, r, *fill);
	if (border)
		DrawCircleLinesV(v, r, *border);
}
