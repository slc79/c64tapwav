#include "interpolate.h"

double lanczos_table[(LANCZOS_RADIUS * 2) * LANCZOS_RESOLUTION];

void make_lanczos_weight_table()
{
	for (int i = 0; i < (LANCZOS_RADIUS * 2) * LANCZOS_RESOLUTION; ++i) {
		float x = double(i) / LANCZOS_RESOLUTION - LANCZOS_RADIUS;
		lanczos_table[i] = lanczos_weight(x);
	}
}
