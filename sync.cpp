// A program that tries to sync up the two channels, using Viterbi decoding
// to find the most likely misalignment as it changes throughout the file.

#define NUM_THEORIES 200
#define THEORY_FROM -20.0  /* in samples */
#define THEORY_TO 20.0  /* in samples */
#define SWITCH_COST 1000.0  /* pretty arbitrary */
#define _CRT_SECURE_NO_WARNINGS

#if defined(_MSC_VER)
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif

#include <stdio.h>
#include <math.h>
#include <vector>
#include <algorithm>

#include "interpolate.h"

#define BUFSIZE 4096

struct wav_file{
	unsigned int chunk_id;
	unsigned int chunk_size;
	unsigned int format;
	unsigned int subchunk1_id;
	unsigned int subchunk1_size;
	unsigned short audio_format;
	unsigned short num_channels;
	unsigned int sample_rate;
	unsigned int byte_rate;
	unsigned short block_align;
	unsigned short bits_per_sample;
	unsigned int subchunk2_id;
	unsigned int subchunk2_size;
	char *data;
};

struct stereo_sample {
	short left, right;
};
struct float_stereo_sample {
	float left, right;
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

struct hypothesis {
	int id;
	double cost;
	hypothesis *prev;
};

int wav_sanitycheck(wav_file wav)
{

	printf("Sample rate : %d Hz\n", wav.sample_rate);
	printf("Bit depth   : %d bit\n", wav.bits_per_sample);
	printf("Channels    : %d\n\n", wav.num_channels);

	if (strncmp((char*)&wav.chunk_id, "RIFF", 4) != 0)
	{
		printf("Not a RIFF file\n");
		return -1;
	}
	
	if (strncmp((char*)&wav.subchunk1_id, "fmt ", 4) != 0)
	{
		printf("Unexpected subchunk\n");
		return -1;
	}
	if (wav.num_channels != 2)
	{
		printf("Only stereo supported\n");
		return -1;
	}

	if (wav.audio_format != 1)
	{
		printf("Only PCM supported\n");
		return -1;
	}

	if (wav.bits_per_sample != 16)
	{
		printf("Only 16 bit audio supported\n");
		return -1;
	}

	if (wav.subchunk1_size != 16)
	{
		printf("Unsupported wav-format %d", wav.subchunk1_size);
		return -1;
	}

	return 0;
}
int main(int argc, char **argv)
{
	if (argc < 3)
	{
		printf("Usage: sync <source.wav> <destination.wav>\n");
		exit(0);
	}

	FILE *ftr;

	make_lanczos_weight_table();
	std::vector<stereo_sample> pcm;

	float sample_rate = 48000.0f;
	unsigned int filesize = 0;

	wav_file source_file;
	wav_file dest_file;
	ftr = fopen(argv[1], "rb");
	fseek(ftr, 0, SEEK_END);
	filesize = ftell(ftr);
	rewind(ftr);

	if(filesize > 44)
		fread(&source_file, sizeof(char), 44, ftr);
	else
	{
		printf("Invalid file\n");
		exit(0);
	}

	if (wav_sanitycheck(source_file) == -1)
	{
		printf("Loading terminated\n");
		exit(0);
	}

	fseek(ftr, 44, SEEK_SET);
	source_file.data = (char*)malloc((sizeof(char)*filesize) - 44);
	fread(source_file.data, sizeof(char), filesize - 44, ftr);
	fclose(ftr);

	unsigned int number_of_samples;
	sample_rate = (float)source_file.sample_rate;

	for (unsigned int i = 0;i < filesize-44;i+=4)
	{
		stereo_sample read_sample;
		read_sample.left = (short)((unsigned char)source_file.data[i + 1]) << 8 | (unsigned char)source_file.data[i];
		read_sample.right = (short)((unsigned char)source_file.data[i + 3] << 8 | (unsigned char)source_file.data[i + 2]);
		pcm.push_back(read_sample);
	}

	free(source_file.data);
	source_file.data = NULL;

	double sum_left = 0.0, sum_right = 0.0;
	for (unsigned i = 0; i < pcm.size(); ++i) {
		sum_left += pcm[i].left;
		sum_right += pcm[i].right;
	}
	double mean_left = sum_left / pcm.size();
	double mean_right = sum_right / pcm.size();

	//fprintf(stderr, "Mean: L=%f R=%f\n", mean_left, mean_right);

	double sum2_left = 0.0, sum2_right = 0.0;
	for (unsigned i = 0; i < pcm.size(); ++i) {
		sum2_left += (pcm[i].left - mean_left) * (pcm[i].left - mean_left);
		sum2_right += (pcm[i].right - mean_right) * (pcm[i].right - mean_right);
	}
	double var_left = sum2_left / (pcm.size() - 1);
	double var_right = sum2_right / (pcm.size() - 1);

	//fprintf(stderr, "Stddev: L=%f R=%f\n", sqrt(var_left), sqrt(var_right));

	double inv_sd_left = 1.0 / sqrt(var_left);
	double inv_sd_right = 1.0 / sqrt(var_right);

	std::vector<float_stereo_sample> norm;
	norm.resize(pcm.size());

	for (unsigned i = 0; i < pcm.size(); ++i) {
		norm[i].left = (pcm[i].left - mean_left) * inv_sd_left;
		norm[i].right = (pcm[i].right - mean_right) * inv_sd_right;
	}

#if 0
	double offset = 0.0;
	double old_diff = 0.0;
	for (unsigned i = 0; i < pcm.size(); ++i) {
		double left = (pcm[i].left - mean_left) * inv_sd_left;
		double right = (interpolate(pcm, i + offset) - mean_right) * inv_sd_right;

		double diff = right - left;
		old_diff = old_diff * 0.9999 + diff * 0.0001;
		offset -= 0.1 * old_diff;

		if (i % 100 == 0) {
			fprintf(stderr, "%7.3f: %7.3f [diff=%8.3f lagged diff=%8.3f]\n", i / sample_rate, offset, diff, old_diff);
		}
		printf("%f %f %f\n", i / sample_rate, left, right);
	}
#endif

	double delays[NUM_THEORIES];
	std::vector<hypothesis *> alloc_hypot;
	hypothesis *bases = new hypothesis[NUM_THEORIES];
	alloc_hypot.push_back(bases);

	for (int h = 0; h < NUM_THEORIES; ++h) {
		delays[h] = THEORY_FROM + h * (THEORY_TO - THEORY_FROM) / (NUM_THEORIES - 1);
		bases[h].id = h;
		bases[h].cost = 0.0;
		bases[h].prev = NULL;
	}

	fprintf(stderr, "Matching blocks... %7.2f", 0.0);
	hypothesis *prev_hyp = bases;
	size_t total_end = pcm.size();
	//size_t total_end = 441000;
	for (unsigned i = 0; i < total_end; i += BUFSIZE) {
		fprintf(stderr, "\b\b\b\b\b\b\b%7.2f", i / sample_rate);
		size_t end = std::min<size_t>(i + BUFSIZE, total_end);
	
		hypothesis *hyp = new hypothesis[NUM_THEORIES];
		alloc_hypot.push_back(hyp);

		// evaluate all hypotheses
		for (int h = 0; h < NUM_THEORIES; ++h) {
			double sum = 0.0;
			double d = delays[h];
			for (unsigned s = i; s < end; ++s) {
				double left = norm[s].left;
				double right = linear_interpolate_right(norm, s + d);
				double diff = (right - left) * (right - left);
				sum += diff;
			}

			double best_cost = HUGE_VAL;
			hypothesis *best_prev = NULL;
			for (int hp = 0; hp < NUM_THEORIES; ++hp) {
				double switch_cost = SWITCH_COST * fabs(delays[hp] - delays[h]);
				double cost = prev_hyp[hp].cost + sum + switch_cost;
				if (best_prev == NULL || cost < best_cost) {
					best_cost = cost;
					best_prev = &prev_hyp[hp];
				}
			}

			hyp[h].id = h;
			hyp[h].cost = best_cost;
			hyp[h].prev = best_prev;
		}

		prev_hyp = hyp;
	}
	fprintf(stderr, "\b\b\b\b\b\b\b%7.2f\n", total_end / sample_rate);

	// best winner
	double best_cost = HUGE_VAL;
	hypothesis *best_hyp = NULL;
	for (int h = 0; h < NUM_THEORIES; ++h) {
		if (best_hyp == NULL || prev_hyp[h].cost < best_cost) {
			best_cost = prev_hyp[h].cost;
			best_hyp = &prev_hyp[h];
		}
	}

	// trace the path backwards
	std::vector<double> best_path;
	while (best_hyp != NULL) {
		best_path.push_back(delays[best_hyp->id]);
		best_hyp = best_hyp->prev;
	}

	reverse(best_path.begin(), best_path.end());

	// save some RAM

	FILE *fp;

	norm = std::vector<float_stereo_sample>();
	for (unsigned i = 0; i < alloc_hypot.size(); ++i) {
		delete[] alloc_hypot[i];
	}

	fprintf(stderr, "Stretching right channel to match left... %7.2f%%", 0.0);
	double inv_sd = sqrt(2.0) / sqrt(var_left + var_right);
	std::vector<stereo_sample> aligned_pcm;
	std::vector<short> mono_pcm;
	aligned_pcm.resize(total_end);
	mono_pcm.resize(total_end);
	for (unsigned i = 0; i < total_end; ++i) {
		double d = lanczos_interpolate(best_path, i / double(BUFSIZE));
		int left = pcm[i].left;
		int right = lanczos_interpolate_right(pcm, i + d);

		aligned_pcm[i].left = left;
		aligned_pcm[i].right = clip(right);

		mono_pcm[i] = clip((int)(inv_sd * 4096.0 * (left + right)));

		if (i % 4096 == 0) {
			fprintf(stderr, "\b\b\b\b\b\b\b\b%7.2f%%", 100.0 * i / total_end);
		}
	}
	fprintf(stderr, "\b\b\b\b\b\b\b%7.2f%%\n", 100.0);

	dest_file.chunk_id = 0x46464952; // RIFF
	dest_file.chunk_size = 36 + (sizeof(short) * mono_pcm.size());
	dest_file.format = 0x45564157; // WAVE
	dest_file.subchunk1_id = 0x20746d66; // fmt
	dest_file.subchunk1_size = 16;
	dest_file.audio_format = 1;
	dest_file.num_channels = 1;
	dest_file.bits_per_sample = 16;
	dest_file.sample_rate = (unsigned int)sample_rate;
	dest_file.byte_rate = (dest_file.sample_rate * dest_file.num_channels) + (dest_file.bits_per_sample / 8);
	dest_file.block_align = dest_file.num_channels * (dest_file.bits_per_sample / 8);
	dest_file.subchunk2_id = 0x61746164;
	dest_file.subchunk2_size = sizeof(short) * mono_pcm.size();
	
	fprintf(stderr, "Writing combined mono track...\n");
	fp = fopen(argv[2], "wb");
	fwrite(&dest_file, sizeof(char) * 44, 1, fp);
	fwrite(mono_pcm.data(), sizeof(short) * mono_pcm.size(), 1, fp);
	fclose(fp);

	fprintf(stderr, "All done.\n");

#if 0

	for (int sec = 0; sec < 2400; ++sec) {
		double from_offset = -128.0;
		double to_offset = 128.0;
		double best_offset = HUGE_VAL;
		for (int iter = 0; iter < 5; ++iter) {	
			//printf("from %f to %f\n", from_offset, to_offset);
			double best_sum = HUGE_VAL;
		
			for (int k = 0; k < 100; ++k) {
				double offset = from_offset + k * 0.01 * (to_offset - from_offset);
				double sum = 0.0;
				for (unsigned i = sec * 4410; i < sec * 4410 + 4410; ++i) {
				//for (unsigned i = 100000; i < 10NUM_THEORIES0; ++i) {
					double left = norm[i].left;
					double right = interpolate(norm, i + offset);
					//double right = (norm[i + offset].right - mean_right) * inv_sd_right;
					double diff = (right - left) * (right - left);
					sum += diff;
				}
		
				if (sum < best_sum) {
					best_sum = sum;
					best_offset = offset;
				}

				//printf("%f %f\n", offset, sqrt(sum / NUM_THEORIES00.0f));
				//fflush(stdout);
			}
		
			double range = 0.1 * (to_offset - from_offset);
			from_offset = best_offset - range;
			to_offset = best_offset + range;
		}
		printf("%d %f\n", sec, best_offset);
		fflush(stdout);
	}
#endif
	
}
