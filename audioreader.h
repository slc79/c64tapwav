#ifndef _AUDIOREADER_H
#define _AUDIOREADER_H 1

#include <vector>

#include <stdint.h>

bool read_audio_file(const char *filename, std::vector<int16_t> *samples);

#endif  // !defined(_AUDIOREADER_H)
