#pragma once
#include "Header.h"
using namespace std::literals::complex_literals;

namespace vis {
	inline Cf downgrade(C c) { return Cf((float)c.real(), (float)c.imag()); }
	inline Vector2 v2(Cf c) { return Vector2{ c.real(), c.imag() }; }
	inline Vector2 v2(C c) { return v2(downgrade(c)); }
	void plot(Cf v, Color c = BLACK);
	void circle(Cf o, float r, Color const* fill, Color const* border);
};
