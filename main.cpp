#include <iostream>
#include <random>
#include <cmath>
#include <time.h>

#define SAMPLE_SIZE (16)
#define ITERATIONS (100000)

std::uniform_real_distribution<double> dist(0.0, 1.0);

double generator()
{
	std::random_device rd; // nondeterministic seed
	std::mt19937 gen(rd());

	return dist(gen);
}

#define PI (3.14159265358979323846264338327950288419716939937510582097494459230781640628)
#define a1 (-3.969683028665376e+01)
#define a2 (2.209460984245205e+02)
#define a3 (-2.759285104469687e+02)
#define a4 (1.383577518672690e+02)
#define a5 (-3.066479806614716e+01)
#define a6 (2.506628277459239e+00)

#define b1 (-5.447609879822406e+01)
#define b2 (1.615858368580409e+02)
#define b3 (-1.556989798598866e+02)
#define b4 (6.680131188771972e+01)
#define b5 (-1.328068155288572e+01)

#define c1 (-7.784894002430293e-03))
#define c2 (-3.223964580411365e-01)
#define c3 (-2.400758277161838e+00)
#define c4 (-2.549732539343734e+00)
#define c5 (4.374664141464968e+00)
#define c6 (2.938163982698783e+00)

#define d1 (7.784695709041462e-03)
#define d2 (3.224671290700398e-01)
#define d3 (2.445134137142996e+00)
#define d4 (3.754408661907416e+00)

#define p_low (0.02425)
#define p_high (0.97575) // 1 - p_low

// https://web.archive.org/web/20150910144400/http://home.online.no/~pjacklam/notes/invnorm
double ncdinv(double p)
{
	double q;
	double r;
	double e;
	double u;
	double x; // out

	if (p < 0)
	{
		return -INFINITY;
	}

	if (p > 1)
	{
		return INFINITY;
	}

	if (p < p_low)
	{
		q = sqrt(-2 * log(p));
		x = ((((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
	} else if (p_low <= p && p <= p_high)
	{
		q = p - 0.5;
		r = q * q;
		x = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /
			(((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1);
	} else if (p > p_high)
	{
		q = sqrt(-2 * log(1 - p));
		x = -((((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
	}

	// Halley rational method for refinement
	//e = Phi(x) - p;								// [use this line
	e = 0.5 * erfc(-x / sqrt(2)) - p;        	// OR this; not both]
	u = e * sqrt(2 * PI) * exp(x * x / 2);
	x = x - u / (1 + x * u / 2);

	return x;
}

void simulate_distribution()
{
	timespec start, end;

	clock_gettime(CLOCK_MONOTONIC, &start);

	unsigned int distribution[SAMPLE_SIZE] = {0};

	for (int i = 0; i < ITERATIONS; ++i)
	{
		int index = round(ncdinv(generator()) + SAMPLE_SIZE * 0.5);

		if (index < 0 || index > SAMPLE_SIZE - 1)
		{
			// could generate outside of sample size (very unlikely
			std::cout << "Hit outside of sample size: " << index << std::endl;
			continue;
		}

		++distribution[index];
	}

	clock_gettime(CLOCK_MONOTONIC, &end);
	std::cout << "Completed in: " << (end.tv_sec - start.tv_sec + (end.tv_nsec - start.tv_nsec) / 1'000'000'000) << " s" << std::endl;

	for (int i = 0; i < SAMPLE_SIZE; ++i)
	{
		std::cout << distribution[i] << " ";
	}

	std::cout.flush();
}

int main()
{
	simulate_distribution();
	return 0;
}
