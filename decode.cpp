#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <vector>
#include <algorithm>

#include "audioreader.h"
#include "interpolate.h"
#include "tap.h"

#define BUFSIZE 4096
#define HYSTERESIS_LIMIT 3000
#define C64_FREQUENCY 985248

#define SYNC_PULSE_START 1000
#define SYNC_PULSE_END 20000
#define SYNC_PULSE_LENGTH 378.0
#define SYNC_TEST_TOLERANCE 1.10

// between [x,x+1]
double find_zerocrossing(const std::vector<float> &pcm, int x)
{
	if (pcm[x] == 0) {
		return x;
	}
	if (pcm[x + 1] == 0) {
		return x + 1;
	}

	assert(pcm[x + 1] < 0);
	assert(pcm[x] > 0);

	double upper = x;
	double lower = x + 1;
	while (lower - upper > 1e-3) {
		double mid = 0.5f * (upper + lower);
		if (lanczos_interpolate(pcm, mid) > 0) {
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

// Calibrate on the first ~25k pulses (skip a few, just to be sure).
double calibrate(const std::vector<pulse> &pulses) {
	if (pulses.size() < SYNC_PULSE_END) {
		fprintf(stderr, "Too few pulses, not calibrating!\n");
		return 1.0;
	}

	int sync_pulse_end = -1;
	double sync_pulse_stddev = -1.0;

	// Compute the standard deviation (to check for uneven speeds).
	// If it suddenly skyrockets, we assume that sync ended earlier
	// than we thought (it should be 25000 cycles), and that we should
	// calibrate on fewer cycles.
	for (int try_end : { 2000, 4000, 5000, 7500, 10000, 15000, SYNC_PULSE_END }) {
		double sum2 = 0.0;
		for (int i = SYNC_PULSE_START; i < try_end; ++i) {
			double cycles = pulses[i].len * C64_FREQUENCY;
			sum2 += (cycles - SYNC_PULSE_LENGTH) * (cycles - SYNC_PULSE_LENGTH);
		}
		double stddev = sqrt(sum2 / (try_end - SYNC_PULSE_START - 1));
		if (sync_pulse_end != -1 && stddev > 5.0 && stddev / sync_pulse_stddev > 1.3) {
			fprintf(stderr, "Stopping at %d sync pulses because standard deviation would be too big (%.2f cycles); shorter-than-usual trailer?\n",
				sync_pulse_end, stddev);
			break;
		}
		sync_pulse_end = try_end;
		sync_pulse_stddev = stddev;
	}
	fprintf(stderr, "Sync pulse length standard deviation: %.2f cycles\n",
		sync_pulse_stddev);

	double sum = 0.0;
	for (int i = SYNC_PULSE_START; i < sync_pulse_end; ++i) {
		sum += pulses[i].len;
	}
	double mean_length = C64_FREQUENCY * sum / (sync_pulse_end - SYNC_PULSE_START);
	double calibration_factor = SYNC_PULSE_LENGTH / mean_length;
	fprintf(stderr, "Calibrated sync pulse length: %.2f -> %.2f (change %+.2f%%)\n",
		mean_length, SYNC_PULSE_LENGTH, 100.0 * (calibration_factor - 1.0));

	// Check for pulses outside +/- 10% (sign of misdetection).
	for (int i = SYNC_PULSE_START; i < sync_pulse_end; ++i) {
		double cycles = pulses[i].len * calibration_factor * C64_FREQUENCY;
		if (cycles < SYNC_PULSE_LENGTH / SYNC_TEST_TOLERANCE || cycles > SYNC_PULSE_LENGTH * SYNC_TEST_TOLERANCE) {
			fprintf(stderr, "Sync cycle with downflank at %.6f was detected at %.0f cycles; misdetect?\n",
				pulses[i].time, cycles);
		}
	}

	return calibration_factor;
}
	
int main(int argc, char **argv)
{
	make_lanczos_weight_table();
	std::vector<float> pcm;
	int sample_rate;
	if (!read_audio_file(argv[1], &pcm, &sample_rate)) {
		exit(1);
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
	double last_downflank = -1;
	for (unsigned i = 0; i < pcm.size(); ++i) {
		int bit = (pcm[i] > 0) ? 1 : 0;
		if (bit == 0 && last_bit == 1) {
			// Check if we ever go up above HYSTERESIS_LIMIT before we dip down again.
			bool true_pulse = false;
			unsigned j;
			int min_level_after = 32767;
			for (j = i; j < pcm.size(); ++j) {
				min_level_after = std::min<int>(min_level_after, pcm[j]);
				if (pcm[j] > 0) break;
				if (pcm[j] < -HYSTERESIS_LIMIT) {
					true_pulse = true;
					break;
				}
			}

			if (!true_pulse) {
#if 0
				fprintf(stderr, "Ignored down-flank at %.6f seconds due to hysteresis (%d < %d).\n",
					double(i) / sample_rate, -min_level_after, HYSTERESIS_LIMIT);
#endif
				i = j;
				continue;
			} 

			// down-flank!
			double t = find_zerocrossing(pcm, i - 1) * (1.0 / sample_rate);
			if (last_downflank > 0) {
				pulse p;
				p.time = t;
				p.len = t - last_downflank;
				pulses.push_back(p);
			}
			last_downflank = t;
		}
		last_bit = bit;
	}

	double calibration_factor = calibrate(pulses);

	FILE *fp = fopen("cycles.plot", "w");
	std::vector<char> tap_data;
	for (unsigned i = 0; i < pulses.size(); ++i) {
		double cycles = pulses[i].len * calibration_factor * C64_FREQUENCY;
		fprintf(fp, "%f %f\n", pulses[i].time, cycles);
		int len = lrintf(cycles / TAP_RESOLUTION);
		if (i > SYNC_PULSE_END && (cycles < 100 || cycles > 800)) {
			fprintf(stderr, "Cycle with downflank at %.6f was detected at %.0f cycles; misdetect?\n",
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
	fclose(fp);

	tap_header hdr;
	memcpy(hdr.identifier, "C64-TAPE-RAW", 12);
	hdr.version = 1;
	hdr.reserved[0] = hdr.reserved[1] = hdr.reserved[2] = 0;
	hdr.data_len = tap_data.size();

	fwrite(&hdr, sizeof(hdr), 1, stdout);
	fwrite(tap_data.data(), tap_data.size(), 1, stdout);

	// Output a debug raw file with pulse detection points.
	fp = fopen("debug.raw", "wb");
	short one = 32767;
	short zero = 0;
	unsigned pulsenum = 0;
	for (unsigned i = 0; i < pcm.size(); ++i) {
		unsigned next_pulse = (pulsenum >= pulses.size()) ? INT_MAX : int(pulses[pulsenum].time * sample_rate);
		if (i >= next_pulse) {
			fwrite(&one, sizeof(one), 1, fp);
			++pulsenum;
		} else {
			fwrite(&zero, sizeof(zero), 1, fp);
		}
	}
	fclose(fp);
}
