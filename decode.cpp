// Copyright Steinar H. Gunderson <sgunderson@bigfoot.com>
// Licensed under the GPL, v2. (See the file COPYING.)

#define NOMINMAX
#define __AVX__

float filter200hz[88] = {-6.13722842764208e-05,-6.13722842764208e-05,-0.000122744568552842,-0.000184116852829262,-0.000245489137105683,-0.000368233705658525,-0.000490978274211366,-0.000613722842764208,-0.00079783969559347,-0.00104332883269915,-0.00128881796980484,-0.00159567939118694,-0.00196391309684546,-0.00233214680250399,-0.00282312507671536,-0.00337547563520314,-0.00398919847796735,-0.00466429360500798,-0.00540076101632503,-0.0061986007119185,-0.00705781269178839,-0.00803976924021112,-0.00902172578863385,-0.0101264269056094,-0.0112925003068614,-0.0124585737081134,-0.0136860193936418,-0.0149134650791702,-0.0162022830489751,-0.0174911010187799,-0.0187799189885848,-0.0200687369583896,-0.0213575549281944,-0.0225236283294464,-0.0236897017306984,-0.024794402847674,-0.0257763593960967,-0.026696943660243,-0.0274947833558365,-0.0281698784828771,-0.0287222290413649,-0.0291518350312999,-0.0293973241684055,1,-0.0295814410212348,-0.0293973241684055,-0.0291518350312999,-0.0287222290413649,-0.0281085061986007,-0.0274334110715601,-0.0266355713759666,-0.0257149871118203,-0.0247330305633976,-0.023628329446422,-0.02246225604517,-0.021296182643918,-0.0200073646741132,-0.0187799189885848,-0.0174297287345035,-0.0161409107646987,-0.0148520927948938,-0.0136246471093654,-0.012397201423837,-0.011231128022585,-0.010065054621333,-0.00902172578863385,-0.0079783969559347,-0.00705781269178839,-0.00613722842764208,-0.00533938873204861,-0.00460292132073156,-0.00392782619369093,-0.00331410335092672,-0.00282312507671536,-0.00233214680250399,-0.00190254081256904,-0.00159567939118694,-0.00128881796980484,-0.000981956548422732,-0.00079783969559347,-0.000613722842764208,-0.000490978274211366,-0.000368233705658525,-0.000245489137105683,-0.000184116852829262,-0.000122744568552842,-6.13722842764208e-05,-6.13722842764208e-05 };

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <io.h>
#include <fcntl.h>
#ifdef __AVX__
#include <immintrin.h>
#endif
#include <vector>
#include <algorithm>

#include "getopt.h"
#include "audioreader.h"
#include "interpolate.h"
#include "level.h"
#include "tap.h"
#include "filter.h"

#define BUFSIZE 4096
#define C64_FREQUENCY 985248
#define SYNC_PULSE_START 1000
#define SYNC_PULSE_END 20000
#define SYNC_PULSE_LENGTH 378.0
#define SYNC_TEST_TOLERANCE 1.10

// SPSA options
#define NUM_FILTER_COEFF 96
#define NUM_SPSA_VALS (NUM_FILTER_COEFF + 2)
#define NUM_ITER 5000
#define A NUM_ITER/10  // approx
#define INITIAL_A 0.005 // A bit of trial and error...
#define INITIAL_C 0.02  // This too.
#define GAMMA 0.166
#define ALPHA 1.0

static float hysteresis_upper_limit = 0.05;
static float hysteresis_lower_limit = -0.05;
static bool do_calibrate = false;
static bool output_cycles_plot = false;
static bool do_crop = false;
static bool do_range = false;
static bool do_invert = false;
static bool do_tap32 = false;
static float crop_start = 0.0f, crop_end = HUGE_VAL;
static bool apply_range = false;
static float range_start = 0.0f, range_end = HUGE_VAL;

static bool use_fir_filter = false;
static float filter_coeff[NUM_FILTER_COEFF] = { 1.0f };  // The rest is filled with 0.
static bool use_rc_filter = false;
static float rc_filter_freq;
static bool output_filtered = false;

static bool quiet = false;
static bool do_auto_level = false;
static bool output_leveled = false;
static std::vector<float> train_snap_points;
static bool do_train = false;

// The frequency to filter on (for do_auto_level), in Hertz.
// Larger values makes the compressor react faster, but if it is too large,
// you'll ruin the waveforms themselves.
static float auto_level_freq = 200.0;

// The minimum estimated sound level (for do_auto_level) at any given point.
// If you decrease this, you'll be able to amplify really silent signals
// by more, but you'll also increase the level of silent (ie. noise-only) segments,
// possibly caused misdetected pulses in these segments.
static float min_level = 0.05f;

// search for the value <limit> between [x,x+1]
template<bool fast>
double find_crossing(const std::vector<float> &pcm, int x, float limit)
{
	if (fast) {
		// Do simple linear interpolation.
		return x + (limit - pcm[x]) / (pcm[x + 1] - pcm[x]);
	} else {
		// Binary search for the zero crossing as given by Lanczos interpolation.
		double upper = x;
		double lower = x + 1;
		while (lower - upper > 1e-3) {
			double mid = 0.5f * (upper + lower);
			if (lanczos_interpolate(pcm, mid) > limit) {
				upper = mid;
			} else {
				lower = mid;
			}
		}

		return 0.5f * (upper + lower);
	}
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
			fprintf(stderr, "Cycle with downflank at %.6f was detected at %.0f cycles\n",
					pulses[i].time, cycles);
		}


		if (len <= 255 && do_tap32 == false) {
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

/*

For internal testing only.

void output_dmp(const std::vector<pulse>& pulses, double calibration_factor)
{
	std::vector<short> tap_data;
	dmp_header hdr;
	memcpy(hdr.identifier, "DC2N-TAP-RAW", 12);
	hdr.version = 0;
	hdr.machine = 0;
	hdr.resolution = 16;
	hdr.frequency = C64_FREQUENCY;
	hdr.videosystem = 0;

	unsigned int pause_count = 0;

	fwrite(&hdr, sizeof(hdr), 1, stdout);
	
	for (unsigned i = 0; i < pulses.size(); ++i)
	{
		double cycles = pulses[i].len * calibration_factor * C64_FREQUENCY;
		int len = lrintf(cycles);
		if (i > SYNC_PULSE_END && (cycles < 100 || cycles > 800)) 
		{
			pause_count++;
			fprintf(stderr, "Cycle with downflank at %.6f was detected at %.0f cycles\n",
				pulses[i].time, cycles);
		}

		if(len < 0xffff)
			tap_data.push_back(len);
		else
		{
			int overflow_count = len / 0xffff;
			for (int i = 0;i < overflow_count;i++)
				tap_data.push_back(0xffff);
			tap_data.push_back(len % 0xffff);
		}

	}
	
	fprintf(stderr, "\n%d long pulses/pauses detected\n", pause_count);
	fwrite(tap_data.data(), tap_data.size()*2, 1, stdout);
}
*/

static struct option long_options[] = {
	{"auto-level",       0,                 0, 'a' },
	{"auto-level-freq",  required_argument, 0, 'b' },
	{"invert",           0,                 0, 'i' },
	{"output-leveled",   0,                 0, 'A' },
	{"min-level",        required_argument, 0, 'm' },
	{"no-calibrate",     0,                 0, 's' },
	{"plot-cycles",      0,                 0, 'p' },
	{"hysteresis-limit", required_argument, 0, 'l' },
	{"filter",           required_argument, 0, 'f' },
	{"rc-filter",        required_argument, 0, 'r' },
	{"highpass-filter",	 0,					0, 'H' },
	{"output-filtered",  0,                 0, 'F' },
	{"crop",             required_argument, 0, 'c' },
	{"range",            required_argument, 0, 'R' },
	{"train",            required_argument, 0, 't' },
	{"quiet",            0,                 0, 'q' },
	{"help",             0,                 0, 'h' },
	{0,                  0,                 0, 0   }
};

void help()
{
	fprintf(stderr, "decode [OPTIONS] AUDIO-FILE > TAP-FILE\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "  -a, --auto-level             automatically adjust amplitude levels throughout the file\n");
	fprintf(stderr, "  -b, --auto-level-freq        minimum frequency in Hertz of corrected level changes (default 200 Hz)\n");
	fprintf(stderr, "  -i, --invert                 invert audio\n");
	fprintf(stderr, "  -A, --output-leveled         output leveled waveform to leveled.raw\n");
	fprintf(stderr, "  -m, --min-level              minimum estimated sound level (0..1) for --auto-level\n");
	fprintf(stderr, "  -s, --calibrate              try to calibrate on sync pulse length\n");
	fprintf(stderr, "  -p, --plot-cycles            output debugging info to cycles.plot\n");
	fprintf(stderr, "  -l, --hysteresis-limit U[:L] change amplitude threshold for ignoring pulses (-1..1)\n");
	fprintf(stderr, "  -f, --filter C1:C2:C3:...    specify FIR filter (up to %d coefficients)\n", NUM_FILTER_COEFF);
	fprintf(stderr, "  -r, --rc-filter FREQ         send signal through a highpass RC filter with given frequency (in Hertz)\n");
	fprintf(stderr, "  -H, --highpass-filter        use a pre-defined 200Hz highpass filter (do not combine with -r or -f)\n");
	fprintf(stderr, "  -F, --output-filtered        output filtered waveform to filtered.raw\n");
	fprintf(stderr, "  -c, --crop START[:END]       use only the given part of the file\n");
	fprintf(stderr, "  -R, --range START[:END]      applies FIR-filter only to specified range\n");
	fprintf(stderr, "  -t, --train LEN1:LEN2:...    train a filter for detecting any of the given number of cycles\n");
	fprintf(stderr, "                               (implies --no-calibrate and --quiet unless overridden)\n");
	fprintf(stderr, "  -T, --tap32                  outputs a TAP file with true pulse lengths\n");
	fprintf(stderr, "  -q, --quiet                  suppress some informational messages\n");
	fprintf(stderr, "  -h, --help                   display this help, then exit\n");
	exit(1);
}

void parse_options(int argc, char **argv)
{
	for ( ;; ) {
		int option_index = 0;
		int c = getopt_long(argc, argv, "ab:iAm:spl:f:Hr:Fc:R:t:qhT", long_options, &option_index);
		if (c == -1)
			break;

		switch (c) {
		case 'a':
			do_auto_level = true;
			break;

		case 'i':
			do_invert = true;
			break;

		case 'b':
			auto_level_freq = atof(optarg);
			break;

		case 'A':
			output_leveled = true;
			break;

		case 'm':
			min_level = atof(optarg);
			break;

		case 's':
			do_calibrate = true;
			break;
		case 'p':
			output_cycles_plot = true;
			break;

		case 'l': {
			const char *hyststr = strtok(optarg, ": ");
			hysteresis_upper_limit = atof(hyststr);
			hyststr = strtok(NULL, ": ");
			if (hyststr == NULL) {
				hysteresis_lower_limit = -hysteresis_upper_limit;
			} else {
				hysteresis_lower_limit = atof(hyststr);
			}
			break;
		}

		case 'f': {
			const char *coeffstr = strtok(optarg, ": ");
			int coeff_index = 0;
			while (coeff_index < NUM_FILTER_COEFF && coeffstr != NULL) {
				filter_coeff[coeff_index++] = atof(coeffstr);
				coeffstr = strtok(NULL, ": ");
			}
			use_fir_filter = true;
			break;
		}

		case 'H': {
			int coeff_index = 0;
			while (coeff_index < 88) {
				filter_coeff[coeff_index] = filter200hz[coeff_index++];
			}
			use_fir_filter = true;
			break;
		}

		case 'r':
			use_rc_filter = true;
			rc_filter_freq = atof(optarg);
			break;

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

		case 'R': {
			const char *rangestr = strtok(optarg, ":");
			range_start = atof(rangestr);
			rangestr = strtok(NULL, ":");
			if (rangestr == NULL) {
				range_end = HUGE_VAL;
			}
			else {
				range_end = atof(rangestr);
			}
			do_range = true;
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

		case 'T':
			do_tap32 = true;
			break;

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

std::vector<float> do_fir_filter(const std::vector<float>& pcm, const float* filter, int sample_rate)
{
	std::vector<float> filtered_pcm;
	filtered_pcm.resize(pcm.size());
	unsigned i = NUM_FILTER_COEFF;
	unsigned int start_sample = lrintf(range_start * sample_rate);
	unsigned int end_sample = lrintf(range_end * sample_rate);

#ifdef __AVX__
	
	if (do_range == true) {
		unsigned avx_end = i + ((pcm.size() - i) & ~7);
		for ( ; i < avx_end; i += 8) {
			__m256 s = _mm256_setzero_ps();
			for (int j = 0; j < NUM_FILTER_COEFF; ++j) {
				if ((i < end_sample) && (i > start_sample)) {
					__m256 f = _mm256_set1_ps(filter[j]);
					__m256 t;
					t = _mm256_mul_ps(f, _mm256_loadu_ps(&pcm[i - j]));
					s = _mm256_add_ps(t, s);
				}
				else
					s = _mm256_loadu_ps(&pcm[i - j]);
			}
			_mm256_storeu_ps(&filtered_pcm[i], s);
		}
	}
	else {
		unsigned avx_end = i + ((pcm.size() - i) & ~7);
		for (; i < avx_end; i += 8)	{
			__m256 s = _mm256_setzero_ps();
			for (int j = 0; j < NUM_FILTER_COEFF; ++j) {
				__m256 f = _mm256_set1_ps(filter[j]);
				__m256 t;
				//			s = _mm256_fmadd_ps(f, _mm256_load_ps(&pcm[i - j]), s);
				t = _mm256_mul_ps(f, _mm256_loadu_ps(&pcm[i - j]));
				s = _mm256_add_ps(t, s);
			}
			_mm256_storeu_ps(&filtered_pcm[i], s);
		}
	}
#endif

	// Do what we couldn't do with AVX (which is everything for non-AVX machines)
	// as scalar code.
	if (do_range == true) {
		for (; i < pcm.size(); ++i)	{
			float s = 0.0f;
			if ((i > start_sample) && (i < end_sample))	{
				for (int j = 0; j < NUM_FILTER_COEFF; ++j)
					s += filter[j] * pcm[i - j];
			}
			else
				s = pcm[i];

			filtered_pcm[i] = s;
		}
	}
	else
	{
		for (; i < pcm.size(); ++i) {
			float s = 0.0f;
			for (int j = 0; j < NUM_FILTER_COEFF; ++j) {
				s += filter[j] * pcm[i - j];
			}
			filtered_pcm[i] = s;
		}
	}

	if (output_filtered) {
		FILE *fp = fopen("filtered.raw", "wb");
		fwrite(filtered_pcm.data(), filtered_pcm.size() * sizeof(filtered_pcm[0]), 1, fp);
		fclose(fp);
	}

	return filtered_pcm;
}

std::vector<float> do_rc_filter(const std::vector<float>& pcm, float freq, int sample_rate)
{
	// This is only a 6 dB/oct filter, which seemingly works better
	// than the Filter class, which is a standard biquad (12 dB/oct).
	// The b/c calculations come from libnyquist (atone.c);
	// I haven't checked, but I suppose they fall out of the bilinear
	// transform of the transfer function H(s) = s/(s + w).
	std::vector<float> filtered_pcm;
	filtered_pcm.resize(pcm.size());
	const float b = 2.0f - cos(2.0 * M_PI * freq / sample_rate);
	const float c = b - sqrt(b * b - 1.0f);
	float prev_in = 0.0f;
	float prev_out = 0.0f;
	for (unsigned i = 0; i < pcm.size(); ++i) {
		float in = pcm[i];
		float out = c * (prev_out + in - prev_in);
		filtered_pcm[i] = out;
		prev_in = in;
		prev_out = out;
	}

	if (output_filtered) {
		FILE *fp = fopen("filtered.raw", "wb");
		fwrite(filtered_pcm.data(), filtered_pcm.size() * sizeof(filtered_pcm[0]), 1, fp);
		fclose(fp);
	}

	return filtered_pcm;
}

template<bool fast>
std::vector<pulse> detect_pulses(const std::vector<float> &pcm, float hysteresis_upper_limit, float hysteresis_lower_limit, int sample_rate)
{
	std::vector<pulse> pulses;

	// Find the flanks.
	enum State { START, ABOVE, BELOW } state = START;
	double last_downflank = -1;
	for (unsigned i = 0; i < pcm.size(); ++i) {
		if (pcm[i] > hysteresis_upper_limit) {
			state = ABOVE;
		}
		else if (pcm[i] < hysteresis_lower_limit) {
			if (state == ABOVE) {
				// down-flank!
				double t = find_crossing<fast>(pcm, i - 1, hysteresis_lower_limit) * (1.0 / sample_rate) + crop_start;
				if (last_downflank > 0) {
					pulse p;
					p.time = t;
					p.len = t - last_downflank;
					pulses.push_back(p);
				}
				last_downflank = t;
			}
			state = BELOW;
		}
//		else
//			state = ABOVE;
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

std::pair<int, double> find_closest_point(double x, const std::vector<float> &points)
{
	int best_point = 0;
	double best_dist = (x - points[0]) * (x - points[0]);
	for (unsigned j = 1; j < train_snap_points.size(); ++j) {
		double dist = (x - points[j]) * (x - points[j]);
		if (dist < best_dist) {
			best_point = j;
			best_dist = dist;
		}
	}
	return std::make_pair(best_point, best_dist);
}

float eval_badness(const std::vector<pulse>& pulses, double calibration_factor)
{
	double sum_badness = 0.0;
	for (unsigned i = 0; i < pulses.size(); ++i) {
		double cycles = pulses[i].len * calibration_factor * C64_FREQUENCY;
		if (cycles > 2000.0) cycles = 2000.0;  // Don't make pauses arbitrarily bad.
		std::pair<int, double> selected_point_and_sq_dist = find_closest_point(cycles, train_snap_points);
		sum_badness += selected_point_and_sq_dist.second;
	}
	return sqrt(sum_badness / (pulses.size() - 1));
}

void find_kmeans(const std::vector<pulse> &pulses, double calibration_factor, const std::vector<float> &initial_centers)
{
	std::vector<float> last_centers = initial_centers;
	std::vector<float> sums;
	std::vector<float> num;
	sums.resize(initial_centers.size());
	num.resize(initial_centers.size());
	for ( ;; ) {
		for (unsigned i = 0; i < initial_centers.size(); ++i) {
			sums[i] = 0.0f;
			num[i] = 0;
		}
		for (unsigned i = 0; i < pulses.size(); ++i) {
			double cycles = pulses[i].len * calibration_factor * C64_FREQUENCY;
			// Ignore heavy outliers, which are almost always long pauses.
			if (cycles > 2000.0) {
				continue;
			}
			std::pair<int, double> selected_point_and_sq_dist = find_closest_point(cycles, last_centers);
			int p = selected_point_and_sq_dist.first;
			sums[p] += cycles;
			++num[p];
		}
		bool any_moved = false;
		for (unsigned i = 0; i < initial_centers.size(); ++i) {
			if (num[i] == 0) {
				fprintf(stderr, "K-means broke down, can't output new reference training points\n");
				return;
			}
			float new_center = sums[i] / num[i];
			if (fabs(new_center - last_centers[i]) > 1e-3) {
				any_moved = true;
			}
			last_centers[i] = new_center;
		}
		if (!any_moved) {
			break;
		}
	}
	fprintf(stderr, "New reference training points:");
	for (unsigned i = 0; i < last_centers.size(); ++i) {
		fprintf(stderr, " %.3f", last_centers[i]);
	}
	fprintf(stderr, "\n");
}

void spsa_train(const std::vector<float> &pcm, int sample_rate)
{
	float vals[NUM_SPSA_VALS] = { hysteresis_upper_limit, hysteresis_lower_limit, 1.0f };  // The rest is filled with 0.

	float start_c = INITIAL_C;
	double best_badness = HUGE_VAL;

	for (int n = 1; n < NUM_ITER; ++n) {
		float a = INITIAL_A * pow(n + A, -ALPHA);
		float c = start_c * pow(n, -GAMMA);

		// find a random perturbation
		float p[NUM_SPSA_VALS];
		float vals1[NUM_SPSA_VALS], vals2[NUM_SPSA_VALS];
		for (int i = 0; i < NUM_SPSA_VALS; ++i) {
			p[i] = (rand() % 2) ? 1.0 : -1.0;
			vals1[i] = std::max(std::min(vals[i] - c * p[i], 1.0f), -1.0f);
			vals2[i] = std::max(std::min(vals[i] + c * p[i], 1.0f), -1.0f);
		}

		std::vector<pulse> pulses1 = detect_pulses<true>(do_fir_filter(pcm, vals1 + 2, sample_rate), vals1[0], vals1[1], sample_rate);
		std::vector<pulse> pulses2 = detect_pulses<true>(do_fir_filter(pcm, vals2 + 2, sample_rate), vals2[0], vals2[1], sample_rate);
		float badness1 = eval_badness(pulses1, 1.0);
		float badness2 = eval_badness(pulses2, 1.0);

		// Find the gradient estimator
		float g[NUM_SPSA_VALS];
		for (int i = 0; i < NUM_SPSA_VALS; ++i) {
			g[i] = (badness2 - badness1) / (2.0 * c * p[i]);
			vals[i] -= a * g[i];
			vals[i] = std::max(std::min(vals[i], 1.0f), -1.0f);
		}
		if (badness2 < badness1) {
			std::swap(badness1, badness2);
			std::swap(vals1, vals2);
			std::swap(pulses1, pulses2);
		}
		if (badness1 < best_badness) {
			fprintf(stderr, "\nNew best filter (badness=%f): -f", badness1);
			for (int i = 0; i < NUM_FILTER_COEFF; ++i) {
				fprintf(stderr, "%.5f:", vals1[i + 2]);
			}
			fprintf(stderr, " -l%f:%f\n", vals1[0], vals1[1]);
			best_badness = badness1;

			find_kmeans(pulses1, 1.0, train_snap_points);

			if (output_cycles_plot) {
				output_cycle_plot(pulses1, 1.0);
			}
		}
		fprintf(stderr, "%d ", n);
		fflush(stderr);
	}
}

std::vector<float> invert(const std::vector<float> &pcm)
{
	std::vector<float> inverted_pcm;
	inverted_pcm.resize(pcm.size());

	for (int i = 0; i < pcm.size(); i++)
	{
		float sample = pcm[i];
		sample = -sample;
		inverted_pcm[i] = sample;
	}

	return inverted_pcm;
}

int main(int argc, char **argv)
{
	if (argc < 2) {
		help();
		exit(0);
	}

	parse_options(argc, argv);
	_setmode(_fileno(stdout), _O_BINARY);

	make_lanczos_weight_table();
	std::vector<float> pcm;
	int sample_rate;
	if (!read_audio_file(argv[optind], &pcm, &sample_rate)) {
		exit(1);
	}

	if (do_invert) {
		pcm = invert(pcm);
	}

	if (do_crop) {
		pcm = crop(pcm, crop_start, crop_end, sample_rate);
	}

	/*
	if (do_sync) {
		pcm = sync();
	}
	*/

	if (use_fir_filter) {
		pcm = do_fir_filter(pcm, filter_coeff, sample_rate);
	}

	if (use_rc_filter) {
		pcm = do_rc_filter(pcm, rc_filter_freq, sample_rate);
	}

	if (do_auto_level) {
		pcm = level_samples(pcm, min_level, auto_level_freq, sample_rate);
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

	std::vector<pulse> pulses = detect_pulses<false>(pcm, hysteresis_upper_limit, hysteresis_lower_limit, sample_rate);

	double calibration_factor = 1.0;
	if (do_calibrate) {
		calibration_factor = calibrate(pulses);
	}

	if (output_cycles_plot) {
		output_cycle_plot(pulses, calibration_factor);
	}

	output_tap(pulses, calibration_factor);
}
