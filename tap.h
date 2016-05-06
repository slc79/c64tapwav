#ifndef _TAP_H
#define _TAP_H

#define TAP_RESOLUTION 8

struct tap_header {
	char identifier[12];
	char version;
	char reserved[3];
	unsigned int data_len;
};

struct dmp_header {
	char identifier[12];
	char version;
	char machine;
	char videosystem;
	char resolution;
	unsigned int frequency;
};

#endif  // !defined(_TAP_H)
