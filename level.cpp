// A small program to try to even out the audio levels within a file
// (essentially a compressor with infinite lookahead).

#include <stdio.h>
#include <math.h>
#include <vector>
#include <algorithm>

#define BUFSIZE 4096
#define WAVE_FREQ 44100.0

// The frequency to filter on, in Hertz. Larger values makes the
// compressor react faster, but if it is too large, you'll
// ruin the waveforms themselves.
#define LPFILTER_FREQ 50.0

// The minimum estimated sound level at any given point.
// If you decrease this, you'll be able to amplify really silent signals
// by more, but you'll also increase the level of silent (ie. noise-only) segments,
// possibly caused misdetected pulses in these segments.
#define MIN_LEVEL 0.05

// A final scalar to get the audio within approximately the right range.
// Increase to _lower_ overall volume.
#define DAMPENING_FACTOR 5.0

// 6dB/oct per round.
#define FILTER_DEPTH 4

struct stereo_sample {
	short left, right;
};

inline short clip(int x)
{
	if (x < -32768) {
		return x;
	} else if (x > 32767) {
		return 32767;
	} else {
		return short(x);
	}
}

static float a1, a2, b0, b1, b2;
static float d0, d1;

static void filter_init(float cutoff_radians)
{
	float resonance = 1.0f / sqrt(2.0f);
	float sn = sin(cutoff_radians), cs = cos(cutoff_radians);
	float alpha = float(sn / (2 * resonance));

	// coefficients for lowpass filter
        float a0 = 1 + alpha;
	b0 = (1 - cs) * 0.5f;
	b1 = 1 - cs;
	b2 = b0;
        a1 = -2 * cs;
        a2 = 1 - alpha;

	b0 /= a0;
	b1 /= a0;
	b2 /= a0;
	a1 /= a0;
	a2 /= a0;

	// reset filter delays
	d0 = d1 = 0.0f;
}

static float filter_update(float in)
{
	float out = b0*in + d0;
	d0 = b1 * in - a1 * out + d1;
	d1 = b2 * in - a2 * out;
	return out;
}

int main(int argc, char **argv)
{
	std::vector<short> pcm;

	while (!feof(stdin)) {
		short buf[BUFSIZE];
		ssize_t ret = fread(buf, sizeof(short), BUFSIZE, stdin);
		if (ret >= 0) {
			pcm.insert(pcm.end(), buf, buf + ret);
		}
	}
	
	// filter forwards, then backwards (perfect phase filtering)
	std::vector<float> filtered_samples, refiltered_samples;
	filtered_samples.resize(pcm.size());
	refiltered_samples.resize(pcm.size());

	filter_init(M_PI * LPFILTER_FREQ / WAVE_FREQ);
	for (unsigned i = 0; i < pcm.size(); ++i) {
		filtered_samples[i] = filter_update(fabs(pcm[i]));
	}
	filter_init(M_PI * LPFILTER_FREQ / WAVE_FREQ);
	for (unsigned i = pcm.size(); i --> 0; ) {
		refiltered_samples[i] = filter_update(filtered_samples[i]);
	}

	for (int i = 1; i < FILTER_DEPTH; ++i) {
		filter_init(M_PI * LPFILTER_FREQ / WAVE_FREQ);
		for (unsigned i = 0; i < pcm.size(); ++i) {
			filtered_samples[i] = filter_update(refiltered_samples[i]);
		}
		filter_init(M_PI * LPFILTER_FREQ / WAVE_FREQ);
		for (unsigned i = pcm.size(); i --> 0; ) {
			refiltered_samples[i] = filter_update(filtered_samples[i]);
		}
	}

	for (unsigned i = 0; i < pcm.size(); ++i) {
		float f = DAMPENING_FACTOR * std::max(refiltered_samples[i] * (1.0 / 32768.0), MIN_LEVEL);
		short s = clip(lrintf(pcm[i] / f));
		fwrite(&s, sizeof(s), 1, stdout);
	}
}
