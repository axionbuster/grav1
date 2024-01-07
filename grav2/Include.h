#pragma once

#include <functional>
#include <complex>
#include <raylib.h>

typedef std::complex<double> C;
typedef std::complex<float> Cf;

/// <summary>
/// The mathematical constant pi
/// (because Raylib defines "PI" as a 32-bit floating point macro.)
/// </summary>
const double PI64 = std::acos(-1);

/// <summary>
/// Decide whether both components of the vector are finite floating-point numbers.
/// </summary>
inline bool finite(C const& c) { return isfinite(c.real()) && isfinite(c.imag()); }

using namespace std::literals::complex_literals;
