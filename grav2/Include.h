#pragma once

#include <raylib.h>
#include <complex>

typedef std::complex<double> C;
typedef std::complex<float> Cf;

/// <summary>
/// The mathematical constant pi
/// (because Raylib defines "PI" as a 32-bit floating point macro.)
/// </summary>
const double PI64 = std::acos(-1);

using namespace std::literals::complex_literals;
