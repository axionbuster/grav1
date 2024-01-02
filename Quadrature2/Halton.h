#pragma once

namespace halton {
	/// <summary>
	/// Generate a Halton sequence (algorithm is due to Wikipedia) of a given base (`b`).
	/// 
	/// To use Halton low-discrepancy sequences to fill up the unit square
	/// in n-space, use successive prime numbers (2, 3, 5, etc.) and
	/// generate x, y, z, etc. coordinates from these individual streams.
	/// </summary>
	class Halton
	{
		int n{ 0 }, d{ 1 }, x{}, y{}, b{};
	public:
		Halton(int b) : b(b) {}
		double next();
	};
};
