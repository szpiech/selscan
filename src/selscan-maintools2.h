// /* selscan -- a program to calculate EHH-based scans for positive selection in genomes
//    Copyright (C) 2014-2024  Zachary A Szpiech
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 3 of the License, or
//    (at your option) any later version.
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//    You should have received a copy of the GNU General Public License
//    along with this program; if not, write to the Free Software Foundation,
//    Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
// */

// #ifndef __SELSCAN_MAINTOOLS_H__
// #define __SELSCAN_MAINTOOLS_H__

// #include <string>
// #include <map>
// #include <cstdio>
// #include "binom.h"
// #include "param_t.h"
// #include "selscan-data.h"
// #include "selscan-pbar.h"
// #include "hamming_t.h"
// #include "selscan-cli.h"

// using namespace std;


// struct work_order_t
// {
//     int queryLoc;
//     int id;

//     string filename;

//     HaplotypeData *hapData;

//     HaplotypeData *hapData1;
//     HaplotypeData *hapData2;

//     MapData *mapData;

//     double (*calc)(map<string, int> &, int, bool);

//     double *ihs;
//     double *ihhDerivedLeft;
//     double *ihhDerivedRight;
//     double *ihhAncestralLeft;
//     double *ihhAncestralRight;
//     double *freq;

//     double *ihh1;
//     double *freq1;

//     double *ihh2;
//     double *freq2;

//     ofstream *flog;
//     ofstream *fout;
//     Bar *bar;

//     param_t *params;
// };

// struct triplet_t
// {
//     double h1;
//     double h12;
//     double h2dh1;
// };

// void calculatePi(HaplotypeData *hapData, MapData *mapData, int winsize, string outFilename);

// triplet_t calculateSoft(map<string, int> &count, int total);

// void query_locus(void *work_order);
// void query_locus_soft(void *order);

// void calc_ihs(void *work_order);
// void calc_nsl(void *work_order);
// void calc_xpihh(void *work_order);
// void calc_soft_ihs(void *order);

// double calcFreq(HaplotypeData *hapData, int locus, bool unphased);
// int queryFound(MapData *mapData, string query);
// void fillColors(int **hapColor, map<string, int> &hapCount,
//                 string *haplotypeList, int hapListLength,
//                 int currentLoc, int &currentColor, bool left);
// bool familyDidSplit(const string &hapStr, const int hapCount,
//                     int **hapColor, const int nhaps, const int colorIndex,
//                     const int previousLoc, string &mostCommonHap);

// double calculateHomozygosity_Wagh(map<string, int> &count, int total, int derivedCount);

// double calculateHomozygosity(map<string, int> &count, int total, bool ALT);

// #endif