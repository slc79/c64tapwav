#ifndef _INTERPOLATE_H
#define _INTERPOLATE_H 1

#include <math.h>

#include <algorithm>
#include <vector>

#define LANCZOS_RADIUS 30

inline double sinc(double x)
{
	if (fabs(x) < 1e-6) {
		return 1.0f - fabs(x);
	} else {
		return sin(x) / x;
	}
}

inline double lanczos_weight(double x)
{
	if (fabs(x) > LANCZOS_RADIUS) {
		return 0.0f;
	}
	return sinc(M_PI * x) * sinc(M_PI * x / LANCZOS_RADIUS);
}

template<class T>
inline double lanczos_interpolate(const std::vector<T> &pcm, double i)
{
	int lower = std::max<int>(ceil(i - LANCZOS_RADIUS), 0);
	int upper = std::min<int>(floor(i + LANCZOS_RADIUS), pcm.size() - 1);
	double sum = 0.0f;

	for (int x = lower; x <= upper; ++x) {
		sum += pcm[x] * lanczos_weight(i - x);
	}
	return sum;
}

template<class T>
inline double lanczos_interpolate_right(const std::vector<T> &pcm, double i)
{
	int lower = std::max<int>(ceil(i - LANCZOS_RADIUS), 0);
	int upper = std::min<int>(floor(i + LANCZOS_RADIUS), pcm.size() - 1);
	double sum = 0.0f;

	for (int x = lower; x <= upper; ++x) {
		sum += pcm[x].right * lanczos_weight(i - x);
	}
	return sum;
}

template<class T>
inline double linear_interpolate(const std::vector<T> &pcm, double i)
{
	int ii = int(i);
	if (ii < 0 || ii >= int(pcm.size() - 1)) {
		return 0.0;
	}
	double frac = i - ii;

	return pcm[ii] + frac * (pcm[ii + 1] - pcm[ii]);
}

template<class T>
inline double linear_interpolate_right(const std::vector<T> &pcm, double i)
{
	int ii = int(i);
	if (ii < 0 || ii >= int(pcm.size() - 1)) {
		return 0.0;
	}
	double frac = i - ii;

	return pcm[ii].right + frac * (pcm[ii + 1].right - pcm[ii].right);
}
	
#endif  // !defined(_INTERPOLATE_H)
