#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <assert.h>
#include <vector>
#include <algorithm>

#define LANCZOS_RADIUS 30
#define BUFSIZE 4096
#define HYSTERESIS_LIMIT 3000
#define SAMPLE_RATE 44100
#define C64_FREQUENCY 985248
#define TAP_RESOLUTION 8

#define SYNC_PULSE_LENGTH 380.0
#define SYNC_TEST_TOLERANCE 1.10

struct tap_header {
	char identifier[12];
	char version;
	char reserved[3];
	unsigned int data_len;
};

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

double interpolate(const std::vector<short> &pcm, double i)
{
	int lower = std::max<int>(ceil(i - LANCZOS_RADIUS), 0);
	int upper = std::min<int>(floor(i + LANCZOS_RADIUS), pcm.size() - 1);
	double sum = 0.0f;

	for (int x = lower; x <= upper; ++x) {
		sum += pcm[x] * weight(i - x);
	}
	return sum;
}
	
// between [x,x+1]
double find_zerocrossing(const std::vector<short> &pcm, int x)
{
	if (pcm[x] == 0) {
		return x;
	}
	if (pcm[x + 1] == 0) {
		return x + 1;
	}

	assert(pcm[x + 1] > 0);
	assert(pcm[x] < 0);

	double lower = x;
	double upper = x + 1;
	while (upper - lower > 1e-6) {
		double mid = 0.5f * (upper + lower);
		if (interpolate(pcm, mid) > 0) {
			upper = mid;
		} else {
			lower = mid;
		}
	}

	return 0.5f * (upper + lower);
}

struct pulse {
	double time;  // in seconds from start
	double len;   // in seconds
};
	
int main(int argc, char **argv)
{
	std::vector<short> pcm;

	while (!feof(stdin)) {
		short buf[BUFSIZE];
		ssize_t ret = fread(buf, 2, BUFSIZE, stdin);
		if (ret >= 0) {
			pcm.insert(pcm.end(), buf, buf + ret);
		}
	}	

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

	std::vector<pulse> pulses;  // in seconds

	// Find the flanks.
	int last_bit = -1;
	double last_upflank = -1;
	int last_max_level = 0;
	for (unsigned i = 0; i < pcm.size(); ++i) {
		int bit = (pcm[i] > 0) ? 1 : 0;
		if (bit == 1 && last_bit == 0 && last_max_level > HYSTERESIS_LIMIT) {
			// up-flank!
			double t = find_zerocrossing(pcm, i - 1) * (1.0 / SAMPLE_RATE);
			if (last_upflank > 0) {
				pulse p;
				p.time = t;
				p.len = t - last_upflank;
				pulses.push_back(p);
			}
			last_upflank = t;
			last_max_level = 0;
		}
		last_max_level = std::max(last_max_level, abs(pcm[i]));
		last_bit = bit;
	}

	// Calibrate on the first ~25k pulses (skip a few, just to be sure).
	double calibration_factor = 1.0f;
	if (pulses.size() < 20000) {
		fprintf(stderr, "Too few pulses, not calibrating!\n");
	} else {
		double sum = 0.0;
		for (int i = 1000; i < 26000; ++i) {
			sum += pulses[i].len;
		}
		double mean_length = C64_FREQUENCY * sum / 25000.0f;
		calibration_factor = SYNC_PULSE_LENGTH / mean_length;
		fprintf(stderr, "Cycle length: %.2f -> 380.0 (change %+.2f%%)\n",
			mean_length, 100.0 * (calibration_factor - 1.0));
	
		// Check for pulses outside +/- 10% (sign of misdetection).
		for (int i = 1000; i < 25000; ++i) {
			double cycles = pulses[i].len * calibration_factor * C64_FREQUENCY;
			if (cycles < SYNC_PULSE_LENGTH / SYNC_TEST_TOLERANCE || cycles > SYNC_PULSE_LENGTH * SYNC_TEST_TOLERANCE) {
				fprintf(stderr, "Sync cycle with upflank at %.6f was detected at %.0f cycles; misdetect?\n",
					pulses[i].time, cycles);
			}
		}
	}

	std::vector<char> tap_data;
	for (unsigned i = 0; i < pulses.size(); ++i) {
		double cycles = pulses[i].len * calibration_factor * C64_FREQUENCY;
		int len = lrintf(cycles / TAP_RESOLUTION);
		if (i > 15000 && (cycles < 100 || cycles > 800)) {
			fprintf(stderr, "Cycle with upflank at %.6f was detected at %.0f cycles; misdetect?\n",
					pulses[i].time, cycles);
		}
		if (len <= 255) {
			tap_data.push_back(len);
		} else {
			int overflow_len = lrintf(cycles);
			tap_data.push_back(0);
			tap_data.push_back(overflow_len & 0xff);
			tap_data.push_back((overflow_len >> 8) & 0xff);
			tap_data.push_back(overflow_len >> 16);
		}
	}

	tap_header hdr;
	memcpy(hdr.identifier, "C64-TAPE-RAW", 12);
	hdr.version = 1;
	hdr.reserved[0] = hdr.reserved[1] = hdr.reserved[2] = 0;
	hdr.data_len = tap_data.size();

	fwrite(&hdr, sizeof(hdr), 1, stdout);
	fwrite(tap_data.data(), tap_data.size(), 1, stdout);
}
