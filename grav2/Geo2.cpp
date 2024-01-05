#include "Geo2.h"

double Halton::next()
{
	x = d - n;
	if (x == 1)
		n = 1, d *= b;
	else
	{
		y = d / b;
		while (x <= y)
			y /= b;
		n = (b + 1) * y - x;
	}
	return n / (double)d;
}
