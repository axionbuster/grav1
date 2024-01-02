#pragma once
#include "Header.h"
using namespace std::literals::complex_literals;

namespace vis {
	inline Cf downgrade(C c)
	{
		return Cf((float)c.real(), (float)c.imag());
	}
	void plot(Cf v, Color c = BLACK);
	void axes(Cf v, Cf re, Cf im, int n, double label = 0);
	void circle(Cf o, float r, Color const* fill, Color const* border);
};
