#include "hamming_t.h"


int hamming_dist_str(string one, string two)
{
	if(one.compare(two) == 0) return 0;

	int diff = 0;

	for(int i = 0; i < one.length(); i++)
	{
		if(one[i] != two[i]) diff++;
	}

	return diff;
}


int hamming_dist_ptr(short *one, short *two, int length)
{
	if(length == 0) return 0;

	int diff = 0;

	for(int i = 0; i < length; i++)
	{
		if(one[i] != two[i]) diff++;
	}

	return diff;
}

int hamming_dist_ptr(char *one, char *two, int length)
{
	if(length == 0) return 0;

	int diff = 0;

	for(int i = 0; i < length; i++)
	{
		if(one[i] != two[i]) diff++;
	}

	return diff;
}