#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <assert.h>

using namespace std;

#define LANCZOS_RADIUS 10
#define RESOLUTION 256
#define WAVE_FREQ 44100
#define C64_FREQ 985248

// Cutoff frequency of low-pass filter (must be max. WAVE_FREQ / 2 or you
// will get aliasing)
#define LPFILTER_FREQ 22050

#define NORMALIZED_LPFILTER_FREQ (float(LPFILTER_FREQ) / float(WAVE_FREQ))
#define LANCZOS_EFFECTIVE_RADIUS (LANCZOS_RADIUS / (NORMALIZED_LPFILTER_FREQ * 2.0))

#define SIZEOF_HALF_TABLE ((int)((LANCZOS_EFFECTIVE_RADIUS * RESOLUTION)))
static double *integrated_window;

double sinc(double x)
{
	if (fabs(x) < 1e-6) {
		return 1.0f - fabs(x);
	} else {
		return sin(x) / x;
	}
}

double weight(double x)
{
	if (fabs(x) > LANCZOS_EFFECTIVE_RADIUS) {
		return 0.0f;
	}
	float t = 2.0 * M_PI * NORMALIZED_LPFILTER_FREQ * x;
	return sinc(t) * sinc(t / LANCZOS_RADIUS);
}

void make_table()
{
	integrated_window = new double[2 * SIZEOF_HALF_TABLE + 1];

	double sum = 0.0;
	for (int i = 0; i <= 2 * SIZEOF_HALF_TABLE; ++i) {
		float t = (i - SIZEOF_HALF_TABLE) / (float)RESOLUTION;
		integrated_window[i] = sum;
		sum += weight(t) * NORMALIZED_LPFILTER_FREQ * 2.0 / (float)RESOLUTION;
		//printf("%f %f %f\n", t, weight(t), sum);
	}
//	exit(0);
}

// integral from -inf to window
double window_one_sided_integral(double to)
{
	double array_pos = to * RESOLUTION + SIZEOF_HALF_TABLE;

	int whole = int(floor(array_pos));
	double frac = array_pos - whole;
	if (whole < 0) {
		return 0.0;
	}
	if (whole >= 2 * SIZEOF_HALF_TABLE) {
		return 1.0;
	}
	return integrated_window[whole] + frac * (integrated_window[whole + 1] - integrated_window[whole]);
}

double window_integral(double from, double to)
{
	return window_one_sided_integral(to) - window_one_sided_integral(from);
}

struct pulse {
	// all values in samples
	double start, end;
};
vector<pulse> pulses;

int main(int argc, char **argv)
{
	make_table();

	FILE *fp = fopen(argv[1], "rb");
	fseek(fp, 14, SEEK_SET);

	int x = 0;

	while (!feof(fp)) {
		int len = getc(fp);
		int cycles;
		if (len == 0) {
			int a = getc(fp);
			int b = getc(fp);
			int c = getc(fp);
			cycles = a | (b << 8) | (c << 16);
		} else {
			cycles = len * 8;
		}
		pulse p;
		p.start = float(x) * WAVE_FREQ / C64_FREQ;
		p.end = (float(x) + cycles * 0.5) * WAVE_FREQ / C64_FREQ;
		pulses.push_back(p);
		x += cycles;
	}

	int len_cycles = x;
	int len_samples = int(ceil(float(len_cycles) * WAVE_FREQ / C64_FREQ));
	
	fprintf(stderr, "%d pulses, total %.2f seconds (%d samples)\n", pulses.size(), len_cycles / float(C64_FREQ), len_samples);

	int pulse_begin = 0;

	for (int i = 0; i < len_samples; ++i) {
//	for (int i = 50000000; i < 51000000; ++i) {
//		double t = i * 0.01;
		double t = i;
		double sample = 0.0;
		for (unsigned j = pulse_begin; j < pulses.size(); ++j) {
			if (t - pulses[j].end > LANCZOS_EFFECTIVE_RADIUS) {
				++pulse_begin;
				continue;
			}
			if (pulses[j].start - t > LANCZOS_EFFECTIVE_RADIUS) {
				break;
			}
			float contribution = window_integral(pulses[j].start - t, pulses[j].end - t);
			sample += contribution;
		}
		printf("%f\n", sample);
//		short s = lrintf((sample - 0.5f) * 16384.0f);
//		fwrite(&s, 2, 1, stdout);
	}
}
