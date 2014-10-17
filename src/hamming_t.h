#ifndef __HAMMING_T_H__
#define __HAMMING_T_H__

#include <string>

using namespace std;

int hamming_dist_str(string one, string two);
int hamming_dist_ptr(short *one, short *two, int length);
int hamming_dist_ptr(char *one, char *two, int length);

#endif