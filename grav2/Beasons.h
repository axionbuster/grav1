#pragma once
#include "Include.h"

/// <summary>
/// Beason's Runge-Kutta integrator with some modifications.
/// 
/// Runge-Kutta schemes integrate first derivatives only, though
/// I have second derivatives of position (accelerations), and need to
/// integrate them.
/// 
/// Will Beason's method adapts an RK (Runge-Kutta) scheme for second-derivative problems.
/// This one is based on his second blog post in which he claims it has the same order
/// of error as the regular RK schemes from which it is adapted
/// (Note: Beason's method is a general way of adapting an existing RK scheme to
/// make them integrate second derivatives, instead of the first derivatives.)
/// 
/// I adapted his method to the Bogacki-Shampine third-order (global error O(h^3),
/// where h is the duration of each time step) integrator, which is an
/// RK scheme. Bogacki-Shampine integration can estimate the error in each
/// step so that I can attempt the integration again and/or adjust the time step.
/// It additionally fully calculates the acceleration needed for the subsequent
/// time step, which means that time is saved because it can just be recycled.
/// 
/// Validation pending...
/// </summary>
namespace beasons
{
	/// <summary>
	/// A subroutine to take the position and velocity, in this order,
	/// and then return the acceleration.
	/// </summary>
	typedef std::function<C(C, C)> ReckonSecondDerivative;

	/// <summary>
	/// The results of an integration step.
	/// 
	/// The "strong" values are the suggested values; "weak" values
	/// can be compared with the respective strong values to
	/// estimate error. This allows the step size to be adapted
	/// dynamically (functionally not included).
	/// 
	/// Additionally, to save some compute, the acceleration for the
	/// next time step is recorded. I can copy this `y2` value
	/// for the next acceleration calculation (the `y2`
	/// parameter in the `step` function). This is a special property
	/// (known as FSAL --- first same as last property)
	/// of the Bogacki-Shampine method.
	/// </summary>
	struct BeasonsResults
	{
		C y0_strong, y0_weak;
		C y1_strong, y1_weak;
		C y2;
	};

	/// <summary>
	/// Evolve both y and the first derivative of y.
	/// </summary>
	/// <param name="h">Step size</param>
	/// <param name="f">How to compute the second derivative of y</param>
	/// <param name="y0">y</param>
	/// <param name="y1">The first derivative of y</param>
	/// <param name="y2">The second derivative of y, or a non-finite-floating-point number if must be calculated here</param>
	/// <returns>Evolved values of y and the first derivative of y.
	/// The return value contains two estimations of each value,
	/// which can be compared with each other to
	/// estimate the error. It also contains the acceleration
	/// for the next time step, which can be directly plugged in
	/// for the next invocation of this function (into the parameter `y2`, of course).</returns>
	BeasonsResults beason_bogacki_shampine(double h, ReckonSecondDerivative const& f, C y0, C y1, C y2 = 1. / 0.);
}
