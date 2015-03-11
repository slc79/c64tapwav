#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <getopt.h>
#include <vector>
#include <algorithm>

#include "audioreader.h"
#include "interpolate.h"
#include "level.h"
#include "tap.h"

#define BUFSIZE 4096
#define C64_FREQUENCY 985248
#define SYNC_PULSE_START 1000
#define SYNC_PULSE_END 20000
#define SYNC_PULSE_LENGTH 378.0
#define SYNC_TEST_TOLERANCE 1.10

// SPSA options
#define NUM_FILTER_COEFF 32
#define NUM_ITER 5000
#define A NUM_ITER/10  // approx
#define INITIAL_A 0.005 // A bit of trial and error...
#define INITIAL_C 0.02  // This too.
#define GAMMA 0.166
#define ALPHA 1.0

static float hysteresis_limit = 3000.0 / 32768.0;
static bool do_calibrate = true;
static bool output_cycles_plot = false;
static bool use_filter = false;
static bool do_crop = false;
static float crop_start = 0.0f, crop_end = HUGE_VAL;
static float filter_coeff[NUM_FILTER_COEFF] = { 1.0f };  // The rest is filled with 0.
static bool output_filtered = false;
static bool quiet = false;
static bool do_auto_level = false;
static bool output_leveled = false;
static std::vector<float> train_snap_points;
static bool do_train = false;

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
	if (!quiet) {
		fprintf(stderr, "Sync pulse length standard deviation: %.2f cycles\n",
			sync_pulse_stddev);
	}

	double sum = 0.0;
	for (int i = SYNC_PULSE_START; i < sync_pulse_end; ++i) {
		sum += pulses[i].len;
	}
	double mean_length = C64_FREQUENCY * sum / (sync_pulse_end - SYNC_PULSE_START);
	double calibration_factor = SYNC_PULSE_LENGTH / mean_length;
	if (!quiet) {
		fprintf(stderr, "Calibrated sync pulse length: %.2f -> %.2f (change %+.2f%%)\n",
			mean_length, SYNC_PULSE_LENGTH, 100.0 * (calibration_factor - 1.0));
	}

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

void output_tap(const std::vector<pulse>& pulses, double calibration_factor)
{
	std::vector<char> tap_data;
	for (unsigned i = 0; i < pulses.size(); ++i) {
		double cycles = pulses[i].len * calibration_factor * C64_FREQUENCY;
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

	tap_header hdr;
	memcpy(hdr.identifier, "C64-TAPE-RAW", 12);
	hdr.version = 1;
	hdr.reserved[0] = hdr.reserved[1] = hdr.reserved[2] = 0;
	hdr.data_len = tap_data.size();

	fwrite(&hdr, sizeof(hdr), 1, stdout);
	fwrite(tap_data.data(), tap_data.size(), 1, stdout);
}

static struct option long_options[] = {
	{"auto-level",       0,                 0, 'a' },
	{"output-leveled",   0,                 0, 'A' },
	{"no-calibrate",     0,                 0, 's' },
	{"plot-cycles",      0,                 0, 'p' },
	{"hysteresis-limit", required_argument, 0, 'l' },
	{"filter",           required_argument, 0, 'f' },
	{"output-filtered",  0,                 0, 'F' },
	{"crop",             required_argument, 0, 'c' },
	{"quiet",            0,                 0, 'q' },
	{"help",             0,                 0, 'h' },
	{0,                  0,                 0, 0   }
};

void help()
{
	fprintf(stderr, "decode [OPTIONS] AUDIO-FILE > TAP-FILE\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "  -a, --auto-level             automatically adjust amplitude levels throughout the file\n");
	fprintf(stderr, "  -A, --output-leveled         output leveled waveform to leveled.raw\n");
	fprintf(stderr, "  -s, --no-calibrate           do not try to calibrate on sync pulse length\n");
	fprintf(stderr, "  -p, --plot-cycles            output debugging info to cycles.plot\n");
	fprintf(stderr, "  -l, --hysteresis-limit VAL   change amplitude threshold for ignoring pulses (0..32768)\n");
	fprintf(stderr, "  -f, --filter C1:C2:C3:...    specify FIR filter (up to %d coefficients)\n", NUM_FILTER_COEFF);
	fprintf(stderr, "  -F, --output-filtered        output filtered waveform to filtered.raw\n");
	fprintf(stderr, "  -c, --crop START[:END]       use only the given part of the file\n");
	fprintf(stderr, "  -t, --train LEN1:LEN2:...    train a filter for detecting any of the given number of cycles\n");
	fprintf(stderr, "                               (implies --no-calibrate and --quiet unless overridden)\n");
	fprintf(stderr, "  -q, --quiet                  suppress some informational messages\n");
	fprintf(stderr, "  -h, --help                   display this help, then exit\n");
	exit(1);
}

void parse_options(int argc, char **argv)
{
	for ( ;; ) {
		int option_index = 0;
		int c = getopt_long(argc, argv, "aAspl:f:Fc:t:qh", long_options, &option_index);
		if (c == -1)
			break;

		switch (c) {
		case 'a':
			do_auto_level = true;
			break;

		case 'A':
			output_leveled = true;
			break;

		case 's':
			do_calibrate = false;
			break;

		case 'p':
			output_cycles_plot = true;
			break;

		case 'l':
			hysteresis_limit = atof(optarg) / 32768.0;
			break;

		case 'f': {
			const char *coeffstr = strtok(optarg, ": ");
			int coeff_index = 0;
			while (coeff_index < NUM_FILTER_COEFF && coeffstr != NULL) {
				filter_coeff[coeff_index++] = atof(coeffstr);
				coeffstr = strtok(NULL, ": ");
			}
			use_filter = true;
			break;
		}

		case 'F':
			output_filtered = true;
			break;

		case 'c': {
			const char *cropstr = strtok(optarg, ":");
			crop_start = atof(cropstr);
			cropstr = strtok(NULL, ":");
			if (cropstr == NULL) {
				crop_end = HUGE_VAL;
			} else {
				crop_end = atof(cropstr);
			}
			do_crop = true;
			break;
		}

		case 't': {
			const char *cyclestr = strtok(optarg, ":");
			while (cyclestr != NULL) {
				train_snap_points.push_back(atof(cyclestr));
				cyclestr = strtok(NULL, ":");
			}
			do_train = true;

			// Set reasonable defaults (can be overridden later on the command line).
			do_calibrate = false;
			quiet = true;
			break;
		}

		case 'q':
			quiet = true;
			break;

		case 'h':
		default:
			help();
			exit(1);
		}
	}
}

std::vector<float> crop(const std::vector<float>& pcm, float crop_start, float crop_end, int sample_rate)
{
	size_t start_sample, end_sample;
	if (crop_start >= 0.0f) {
		start_sample = std::min<size_t>(lrintf(crop_start * sample_rate), pcm.size());
	}
	if (crop_end >= 0.0f) {
		end_sample = std::min<size_t>(lrintf(crop_end * sample_rate), pcm.size());
	}
	return std::vector<float>(pcm.begin() + start_sample, pcm.begin() + end_sample);
}

// TODO: Support AVX here.
std::vector<float> do_filter(const std::vector<float>& pcm, const float* filter)
{
	std::vector<float> filtered_pcm;
	filtered_pcm.reserve(pcm.size());
	for (unsigned i = NUM_FILTER_COEFF; i < pcm.size(); ++i) {
		float s = 0.0f;
		for (int j = 0; j < NUM_FILTER_COEFF; ++j) {
			s += filter[j] * pcm[i - j];
		}
		filtered_pcm.push_back(s);
	}

	if (output_filtered) {
		FILE *fp = fopen("filtered.raw", "wb");
		fwrite(filtered_pcm.data(), filtered_pcm.size() * sizeof(filtered_pcm[0]), 1, fp);
		fclose(fp);
	}

	return filtered_pcm;
}

std::vector<pulse> detect_pulses(const std::vector<float> &pcm, int sample_rate)
{
	std::vector<pulse> pulses;

	// Find the flanks.
	int last_bit = -1;
	double last_downflank = -1;
	for (unsigned i = 0; i < pcm.size(); ++i) {
		int bit = (pcm[i] > 0) ? 1 : 0;
		if (bit == 0 && last_bit == 1) {
			// Check if we ever go up above <hysteresis_limit> before we dip down again.
			bool true_pulse = false;
			unsigned j;
			int min_level_after = 32767;
			for (j = i; j < pcm.size(); ++j) {
				min_level_after = std::min<int>(min_level_after, pcm[j]);
				if (pcm[j] > 0) break;
				if (pcm[j] < -hysteresis_limit) {
					true_pulse = true;
					break;
				}
			}

			if (!true_pulse) {
#if 0
				fprintf(stderr, "Ignored down-flank at %.6f seconds due to hysteresis (%d < %d).\n",
					double(i) / sample_rate, -min_level_after, hysteresis_limit);
#endif
				i = j;
				continue;
			} 

			// down-flank!
			double t = find_zerocrossing(pcm, i - 1) * (1.0 / sample_rate) + crop_start;
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
	return pulses;
}

void output_cycle_plot(const std::vector<pulse> &pulses, double calibration_factor)
{
	FILE *fp = fopen("cycles.plot", "w");
	for (unsigned i = 0; i < pulses.size(); ++i) {
		double cycles = pulses[i].len * calibration_factor * C64_FREQUENCY;
		fprintf(fp, "%f %f\n", pulses[i].time, cycles);
	}
	fclose(fp);
}

float eval_badness(const std::vector<pulse>& pulses, double calibration_factor)
{
	double sum_badness = 0.0;
	for (unsigned i = 0; i < pulses.size(); ++i) {
		double cycles = pulses[i].len * calibration_factor * C64_FREQUENCY;
		if (cycles > 2000.0) cycles = 2000.0;  // Don't make pauses arbitrarily bad.
		double badness = (cycles - train_snap_points[0]) * (cycles - train_snap_points[0]);
		for (unsigned j = 1; j < train_snap_points.size(); ++j) {
			badness = std::min(badness, (cycles - train_snap_points[j]) * (cycles - train_snap_points[j]));
		}
		sum_badness += badness;
	}
	return sqrt(sum_badness / (pulses.size() - 1));
}

void spsa_train(std::vector<float> &pcm, int sample_rate)
{
	// Train!
	float filter[NUM_FILTER_COEFF] = { 1.0f };  // The rest is filled with 0.

	float start_c = INITIAL_C;
	double best_badness = HUGE_VAL;

	for (int n = 1; n < NUM_ITER; ++n) {
		float a = INITIAL_A * pow(n + A, -ALPHA);
		float c = start_c * pow(n, -GAMMA);

		// find a random perturbation
		float p[NUM_FILTER_COEFF];
		float filter1[NUM_FILTER_COEFF], filter2[NUM_FILTER_COEFF];
		for (int i = 0; i < NUM_FILTER_COEFF; ++i) {
			p[i] = (rand() % 2) ? 1.0 : -1.0;
			filter1[i] = std::max(std::min(filter[i] - c * p[i], 1.0f), -1.0f);
			filter2[i] = std::max(std::min(filter[i] + c * p[i], 1.0f), -1.0f);
		}

		std::vector<pulse> pulses1 = detect_pulses(do_filter(pcm, filter1), sample_rate);
		std::vector<pulse> pulses2 = detect_pulses(do_filter(pcm, filter2), sample_rate);
		float badness1 = eval_badness(pulses1, 1.0);
		float badness2 = eval_badness(pulses2, 1.0);

		// Find the gradient estimator
		float g[NUM_FILTER_COEFF];
		for (int i = 0; i < NUM_FILTER_COEFF; ++i) {
			g[i] = (badness2 - badness1) / (2.0 * c * p[i]);
			filter[i] -= a * g[i];
			filter[i] = std::max(std::min(filter[i], 1.0f), -1.0f);
		}
		if (badness2 < badness1) {
			std::swap(badness1, badness2);
			std::swap(filter1, filter2);
			std::swap(pulses1, pulses2);
		}
		if (badness1 < best_badness) {
			printf("\nNew best filter (badness=%f):", badness1);
			for (int i = 0; i < NUM_FILTER_COEFF; ++i) {
				printf(" %.5f", filter1[i]);
			}
			best_badness = badness1;
			printf("\n");

			if (output_cycles_plot) {
				output_cycle_plot(pulses1, 1.0);
			}
		}
		printf("%d ", n);
		fflush(stdout);
	}
}

int main(int argc, char **argv)
{
	parse_options(argc, argv);

	make_lanczos_weight_table();
	std::vector<float> pcm;
	int sample_rate;
	if (!read_audio_file(argv[optind], &pcm, &sample_rate)) {
		exit(1);
	}

	if (do_crop) {
		pcm = crop(pcm, crop_start, crop_end, sample_rate);
	}

	if (use_filter) {
		pcm = do_filter(pcm, filter_coeff);
	}

	if (do_auto_level) {
		pcm = level_samples(pcm, sample_rate);
		if (output_leveled) {
			FILE *fp = fopen("leveled.raw", "wb");
			fwrite(pcm.data(), pcm.size() * sizeof(pcm[0]), 1, fp);
			fclose(fp);
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

	if (do_train) {
		spsa_train(pcm, sample_rate);
		exit(0);
	}

	std::vector<pulse> pulses = detect_pulses(pcm, sample_rate);

	double calibration_factor = 1.0;
	if (do_calibrate) {
		calibration_factor = calibrate(pulses);
	}

	if (output_cycles_plot) {
		output_cycle_plot(pulses, calibration_factor);
	}

	output_tap(pulses, calibration_factor);
}
