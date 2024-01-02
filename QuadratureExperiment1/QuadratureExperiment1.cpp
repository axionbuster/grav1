#include <iostream>
#include <fstream>
#include <complex>
#include <functional>
#include <random>
#include <vector>

constexpr int n_monte = 1000;

struct Stat
{
	/// <summary>
	/// The sample average and the second moment.
	/// </summary>
	double avg{}, m2{};
	/// <summary>
	/// Number of items (integer).
	/// </summary>
	double n{};
	/// <summary>
	/// Commit a number.
	/// </summary>
	/// <param name="x"></param>
	void put(double x);
	/// <summary>
	/// Where n &gt; 0, compute the sample standard deviation.
	/// </summary>
	double stdev_sam() { return sqrt(m2 / (n - 1)); }
};

struct Log
{
	/// <summary>
	/// Row.
	/// </summary>
	struct Row
	{
		double x{}, y{};
		/// <summary>
		/// If the first bit is set, in the right circle.
		/// 
		/// If the second bit is set, in the left circle.
		/// </summary>
		int flags{};
		bool right() const { return !!(flags & 1); }
		bool left() const { return !!(flags & 2); }
		bool in() const { return (flags & 3) == 2; }
	};
	std::vector<Row> rows;
	bool closed{};
	/// <summary>
	/// Record a computation.
	/// </summary>
	/// <returns></returns>
	void put(double x, double y, bool in_left, bool in_right)
	{
		if (closed) return;
		Row i{ x, y, (!!in_left << 1) | !!in_right };
		rows.push_back(i);
	}
	/// <summary>
	/// Stop admitting more rows.
	/// </summary>
	void close() { closed = true; }

	/// <summary>
	/// Generate a random integer between 0 (inclusive)
	/// to a given exclusive upper bound.
	/// </summary>
	typedef std::function<int(int)> RandomIndex;

	/// <summary>
	/// Take `n` sample rows.
	/// </summary>
	Log samp_n_repl(int n, RandomIndex ri);
};

struct Ctx
{
	std::mt19937 random;
	Log log;
};

/// <summary>
/// By Monte Carlo quadrature,
/// find the relative area of the lune on the "left circle"
/// (unit circle at the origin) when it intersects with the
/// "right circle" (a circle whose center is at (center, 0)
/// with the radius of `radius`).
/// </summary>
/// <param name="center"></param>
/// <param name="radius"></param>
/// <param name="ctx></param>
/// <returns></returns>
double lune_left(double center, double radius, Ctx& ctx);

static std::mt19937 mkrand()
{
	std::random_device dev;
	return std::mt19937(dev());
}

int main()
{
	std::mt19937 rng = mkrand(), rng2 = mkrand();
	Log log{};
	Ctx ctx{ rng, log };
	constexpr int m = 1000;
	Stat s{};
	for (int i = m - 1; i >= 0; i--)
	{
		s.put(lune_left(1 / sqrt(2), 1 / sqrt(2), ctx));
		// Stop recording after first run.
		ctx.log.close();
	}
	std::cout << "n=" << s.n
		<< ", avg=" << s.avg
		<< ", stdev(sam)=" << s.stdev_sam()
		<< "\n\n";

#define PRINTCSV(f,i,r) f << i << "," << (r).x << "," << (r).y \
<< "," << (r).left() << "," << (r).right() << "," << (r).in() \
<< "\n"

	// Output (full) as CSV
	{
		std::ofstream f;
		f.open("out.csv");
		f << "i,x,y,left?,right?,in?\n";
		for (int i = 1; i <= ctx.log.rows.size(); i++)
		{
			auto const& r = ctx.log.rows[i - 1];
			PRINTCSV(f, i, r);
		}
		f << std::endl;
	}

	// Do the same, but abbreviated.
	{
		constexpr int n_abbrev = 30;
		auto abbrev = ctx.log.samp_n_repl(n_abbrev, [&rng2](int upper)
			{
				std::uniform_int_distribution<> i(0, upper);
				return i(rng2);
			});
		std::ofstream f;
		f.open("outabbr.csv");
		f << "i,x,y,left?,right?,in?\n";
		for (int i = 1; i <= abbrev.rows.size(); i++)
		{
			auto const& r = abbrev.rows[i - 1];
			PRINTCSV(f, i, r);
		}
		f << std::endl;
	}

#undef PRINTCSV

	std::cout << "bye" << std::endl;
}

double lune_left(double c, double r, Ctx& ctx)
{
	using std::min;
	using std::max;
	// Use right-handed coordinates (orientation shouldn't matter).
	double left = min(-1.0, c - r),
		right = max(1.0, c + r),
		top = max(1.0, r),
		bottom = min(-1.0, -r),
		width = right - left,
		height = top - bottom;
	std::uniform_real_distribution<>
		width_distr(left, right),
		height_distr(bottom, top);
	auto X = [&width_distr, &ctx]() { return width_distr(ctx.random); };
	auto Y = [&height_distr, &ctx]() { return height_distr(ctx.random); };
	auto Left = [=](double x, double y) { return x * x + y * y <= 1.0; };
	auto Right = [=](double x, double y)
		{
			double x1 = x - c;
			return x1 * x1 + y * y <= r * r;
		};
	int freq{};
	for (int i = n_monte - 1; i >= 0; i--)
	{
		double x = X(), y = Y();
		bool left = Left(x, y), right = Right(x, y);
		ctx.log.put(x, y, left, right);
		if (left && !right)
			freq++;
	}
	return freq / (double)n_monte;
}

void Stat::put(double x)
{
	// Welford's algorithm.
	double n0 = n++,
		avg0 = avg,
		avg1 = avg += (x - avg) / n;
	m2 += (x - avg0) * (x - avg1);
}

Log Log::samp_n_repl(int n, RandomIndex ri)
{
	std::vector<Log::Row> r;
	while (n-- > 0)
	{
		int k = ri((int)rows.size());
		r.push_back(rows[k]);
	}
	return Log{ r };
}
