#pragma once

// Computation of the area of one side of a lune.

#include <deque>
#include "Header.h"
#include "Halton.h"

namespace lune {
	/// <summary>
	/// Interactive quadrature of the area of the one-sided lune (crescent).
	/// 
	/// When a circle with a given center on the x-axis and radius
	/// intersects with the unit circle at the origin, find the area
	/// of the circle on the right except that part which belongs to
	/// the central unit circle, if that intersecting region exists (or else, zero).
	/// </summary>
	struct Lune
	{
		/// <summary>
		/// Sequence of points that have been recently sampled (FIFO).
		/// 
		/// Back = recent, front = later.
		/// </summary>
		std::deque<C> log;
		/// <summary>
		/// x-coordinate of the center of the right circle;
		/// squared radius of the right circle.
		/// 
		/// The left circle is always the unit circle (radius 1) at the origin.
		/// 
		/// c &gt; 0.
		/// </summary>
		double c{}, rsq{};
		/// <summary>
		/// The number of points satisfying the criteria;
		/// number of points to keep in the log.
		/// </summary>
		int freq{}, cap{};
		/// <summary>
		/// Set up the computation of the area of the right-side of the lune
		/// (crescent) by the intersection of the unit circle at the origin
		/// and a second circle centered at (c, 0) with r > 0 as the radius.
		/// </summary>
		/// <param name="c">x-coordinate of the center of the second circle.</param>
		/// <param name="r">Positive radius of the second circle.</param>
		/// <param name="cap">Number of points to hold in history queue (`log`).</param>
		Lune(double c, double r, int cap = 100);
		/// <summary>
		/// Sample a point and update internal statistics.
		/// </summary>
		void advance();
		/// <summary>
		/// Once at least one point has been sampled,
		/// compute the quadrature by inspection of the internal statistics.
		/// 
		/// (Simple arithmetic).
		/// </summary>
		/// <returns></returns>
		double quadrature() const { return dim * dim * freq / log.size(); }
		/// <summary>
		/// Recall the side length of the bounding square.
		/// </summary>
		/// <returns></returns>
		double dimension() const { return dim; }
		/// <summary>
		/// Recall the x-coordinate of the center of the bounding square.
		/// </summary>
		/// <returns></returns>
		double xmidpoint() const { return m_xmidpoint; }
		/// <summary>
		/// Decide whether the point belongs to the left circle.
		/// </summary>
		/// <param name="p"></param>
		/// <returns></returns>
		static bool left_static(C p) { return std::norm(p) < 1.0; }
		/// <summary>
		/// Decide whether the point belongs to the left circle.
		/// </summary>
		/// <param name="p"></param>
		/// <returns></returns>
		bool left(C p) const { return left_static(p); }
		/// <summary>
		/// Decide whether the point belongs to the right circle.
		/// </summary>
		/// <param name="p"></param>
		/// <returns></returns>
		bool right(C p) const { return std::norm(p - c) < rsq; }
		/// <summary>
		/// Decide whether the point belongs to the right circle
		/// but not the left circle (used for quadrature).
		/// </summary>
		/// <param name="p"></param>
		/// <returns></returns>
		bool in(C p) const { return !left(p) && right(p); }
	private:
		/// <summary>
		/// Internal low-discrepancy sequences used to evenly generate
		/// points in the constructed bounding square, of successive
		/// prime numbers as the "bases" (peculiar details of the algorithm).
		/// 
		/// For each "random" point generated with these:
		/// First a term at h2 is chosen (0 &lt; [h2] &lt; 1), and is turned into the
		/// x-coordinate of that point. The same is done for h3, which becomes
		/// the y-coordinate of that point.
		/// </summary>
		halton::Halton h2, h3;
		/// <summary>
		/// The positive side length of the bounding square.
		/// </summary>
		double dim{};
		/// <summary>
		/// Center of the bounding square (x-coordinate).
		/// </summary>
		double m_xmidpoint{};
	};

	/// <summary>
	/// Computation of the area of a lune of two circles in general position.
	/// </summary>
	struct GenericLune
	{
		/// <summary>
		/// To transform internal coordinates to external coordinates, apply
		/// first the homothety (rotation and scaling), and then the translation, to the internal coordinates.
		/// </summary>
		C homt{ std::polar(1.0, 0.0) }, tr;
		/// <summary>
		/// Internal data structure (intended for public manipulation).
		/// </summary>
		Lune lune;
		/// <summary>
		/// Original circular centers.
		/// </summary>
		C c0, c1;
		/// <summary>
		/// Original radii.
		/// </summary>
		double r0{}, r1{};
		/// <summary>
		/// Prepare the computation in the case of two circles in general position.
		/// </summary>
		/// <param name="c0">Center of the left circle.</param>
		/// <param name="r0">Radius of the left circle (non-zero).</param>
		/// <param name="c1">Center of the right circle.</param>
		/// <param name="r1">Radius of the right circle (non-zero).</param>
		/// <param name="cap">See Lune.</param>
		GenericLune(C c0, double r0, C c1, double r1, int cap = 100);
		/// <summary>
		/// Four corners of the bounding square after the inverse transformation.
		/// </summary>
		struct BoundingSquare
		{
			/// <summary>
			/// Corners 1, 2, 3, and 4, in some contiguous order.
			/// </summary>
			C c[4];
			/// <summary>
			/// Apply the homothety and translation to all four corners.
			/// </summary>
			/// <param name="homt"></param>
			/// <param name="tr"></param>
			void transform(C homt, C tr = 0);
		};
		/// <summary>
		/// Compute the four corners of the bounding square.
		/// </summary>
		/// <returns></returns>
		BoundingSquare bounding() const;
		/// <summary>
		/// Apply the homothety and translation to compute the
		/// original-system coordinates of the point that
		/// the internal `lune` data structure understands.
		/// </summary>
		/// <param name="point">A point understood by `lune`.</param>
		/// <returns>The same point in the original coordinate system.</returns>
		C invert(C point) const;
	};
}
