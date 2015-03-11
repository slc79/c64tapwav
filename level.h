#ifndef _LEVEL_H
#define _LEVEL_H 1

#include <vector>

std::vector<float> level_samples(const std::vector<float> &pcm, float min_level, int sample_rate);

#endif  // !defined(_LEVEL_H)
