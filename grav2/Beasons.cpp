#include "Beasons.h"

using namespace beasons;

/// <summary>
/// Compute the dot product between scalars on the left and
/// vectors on the right, specialized for this integration method.
/// </summary>
/// <param name="left">Vector of coefficients</param>
/// <param name="right">Vector of ... vectors</param>
/// <returns>The dot product</returns>
static C dot(double const left[4], C const right[4])
{
	C z;
	for (int i = 0; i < 4; i++) z += left[i] * right[i];
	return z;
}

/// <summary>
/// Coefficients for the "k" values (here referred to as y2
/// for the second derivative of y).
/// </summary>
static constexpr double A[4][4] = {
	{0, 0, 0, 0},
	{1. / 2, 0, 0, 0},
	{0, 3. / 4, 0, 0},
	{2. / 9, 3. / 9, 4. / 9, 0},
};

/// <summary>
/// The "weak" final coefficient vector. Used to compute error.
/// 
/// Indices are steps (0-indexed).
/// </summary>
static constexpr double bweak[4] = { 2. / 9, 3. / 9, 4. / 9, 0 };

/// <summary>
/// The "strong" final coefficient vector. Suggested for the final value.
/// 
/// Indices are steps (0-indexed).
/// </summary>
static constexpr double bstrong[4] = { 7. / 24, 1. / 4, 1. / 3, 1. / 8 };

/// <summary>
/// Time coefficient vector.
/// 
/// Indices are steps (0-indexed).
/// </summary>
static constexpr double c[4] = { 0, 1. / 2, 3. / 4, 1 };

BeasonsResults
beasons::beason_bogacki_shampine(double h, ReckonSecondDerivative const& f, C y0, C y1, C y2)
{
	// If the second derivative of y is not given, then compute it.
	if (!finite(y2)) y2 = f(y0, y1);

	// Step 0: inputs.
	// Steps 1-3, inclusive: actual work.

	C y0s[4] = { y0, 0, 0, 0 }; // Values of y through the steps.
	C y1s[4] = { y1, 0, 0, 0 }; // First derivatives through the steps.
	C y2s[4] = { y2, 0, 0, 0 }; // Similarly, second derivatives.

	for (int i = 1; i <= 3; i++)
	{
		// Order among the three statements matters.
		// First, zeroth, and then second derivative.

		y1s[i] = y1s[i - 1] + h * dot(A[i], y2s);

		// last term: use y2 at beginning of "step" function call.
		y0s[i] = y0s[i - 1]
			+ 1. / 6 * h
			* (4. * y1s[i - 1] + 2. * y1s[i] + h * c[i] * y2);

		y2s[i] = f(y0s[i], y1s[i]);
	}

	BeasonsResults r;

	r.y1_strong = y1 + dot(bstrong, y2s) * h;
	r.y1_weak = y1 + dot(bweak, y2s) * h;
	r.y0_strong = y0 + dot(bstrong, y1s) * h;
	r.y0_weak = y0 + dot(bweak, y1s) * h;
	r.y2 = y2s[3];

	return r;
}
