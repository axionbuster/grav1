#pragma once

/// <summary>
/// Generation of the Halton low-discrepancy sequence for quasi-Monte Carlo quadrature.
/// </summary>
namespace halton {
	/// <summary>
	/// Generate a Halton sequence (algorithm is due to Wikipedia) of a given base (`b`).
	/// 
	/// To use Halton low-discrepancy sequences to fill up the unit square
	/// in n-space, use successive prime numbers (2, 3, 5, etc.) and
	/// generate x, y, z, etc. coordinates from these individual streams.
	/// </summary>
	struct Halton
	{
		unsigned n{ 0 }, d{ 1 }, x{}, y{}, b{};
	public:
		/// <summary>
		/// Initialize a sequence with the given prime base.
		/// </summary>
		/// <param name="b">Base (prime number).</param>
		Halton(unsigned b) : b(b) {}
		/// <summary>
		/// Extract a term and advance the internal state.
		/// </summary>
		/// <returns></returns>
		double next();
	};
};
