#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <assert.h>
#include <algorithm>

#define LANCZOS_RADIUS 30
#define LEN 813440
#define HYSTERESIS_LIMIT 1000

double sinc(double x)
{
	if (fabs(x) < 1e-6) {
		return 1.0f - fabs(x);
	} else {
		return sin(x) / x;
	}
}

#if 1
double weight(double x)
{
	if (fabs(x) > LANCZOS_RADIUS) {
		return 0.0f;
	}
	return sinc(M_PI * x) * sinc(M_PI * x / LANCZOS_RADIUS);
}
#else
double weight(double x)
{
	if (fabs(x) > 1.0f) {
		return 0.0f;
	}
	return 1.0f - fabs(x);
}
#endif

double interpolate(const short *in, double i)
{
	int lower = std::max(int(ceil(i - LANCZOS_RADIUS)), 0);
	int upper = std::min(int(floor(i + LANCZOS_RADIUS)), LEN - 1);
	double sum = 0.0f;

	for (int x = lower; x <= upper; ++x) {
		sum += in[x] * weight(i - x);
	}
	return sum;
}
	
short in[LEN];

// between [x,x+1]
double find_zerocrossing(int x)
{
	if (in[x] == 0) {
		return x;
	}
	if (in[x + 1] == 0) {
		return x + 1;
	}

	assert(in[x + 1] > 0);
	assert(in[x] < 0);

	double lower = x;
	double upper = x + 1;
	while (upper - lower > 1e-6) {
		double mid = 0.5f * (upper + lower);
		if (interpolate(in, mid) > 0) {
			upper = mid;
		} else {
			lower = mid;
		}
	}

	return 0.5f * (upper + lower);
}

int main(int argc, char **argv)
{
	fread(in, LEN*2, 1, stdin);

#if 0
	for (int i = 0; i < LEN; ++i) {
		in[i] += rand() % 10000;
	}
#endif

#if 0
	for (int i = 0; i < LEN; ++i) {
		printf("%d\n", in[i]);
	}
#endif
	int last_bit = -1;
	double last_upflank = -1;
	int last_max_level = 0;
	for (int i = 0; i < LEN; ++i) {
		int bit = (in[i] > 0) ? 1 : 0;
		if (bit == 1 && last_bit == 0 && last_max_level > HYSTERESIS_LIMIT) {
			// up-flank!
			double t = find_zerocrossing(i - 1) * (123156.0/44100.0);
			if (last_upflank > 0) {
//				fprintf(stderr, "length: %f (0x%x)\n", t - last_upflank, lrintf(t - last_upflank));
				int len = lrintf(t - last_upflank);
				printf("0x%x\n", len);
			}
			last_upflank = t;
			last_max_level = 0;
		}
		last_max_level = std::max(last_max_level, abs(in[i]));
		last_bit = bit;
	}
}
