/* selscan -- a program to calculate EHH-based scans for positive selection in genomes
   Copyright (C) 2014  Zachary A Szpiech

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
*/

#ifndef __XP_IHH_DATA_H__
#define __XP_IHH_DATA_H__
#include <string>
#include <iostream>
#include <fstream>
#include "gzstream.h"

using namespace std;

const int MISSING = -9999;
const char MISSING_CHAR = '9';

struct HaplotypeData
{
    char **data;
    int nhaps;
    int nloci;
};

struct MapData
{
    int *physicalPos;
    double *geneticPos;
    string *locusName;
    int nloci;
    string chr;
};

//allocates the arrays and populates them with -9 or "--" depending on type
MapData *initMapData(int nloci);
void releaseMapData(MapData *data);

//reads in map data and also does basic checks on integrity of format
//returns a populated MapData structure if successful
//throws an exception otherwise
MapData *readMapData(string filename, int expected_loci, bool USE_PMAP);
MapData *readMapDataTPED(string filename, int expected_loci, int expected_haps, bool USE_PMAP);
MapData *readMapDataVCF(string filename, int expected_loci); //Physical positions only

//allocates the 2-d array and populated it with -9
HaplotypeData *initHaplotypeData(unsigned int nhaps, unsigned int nloci);
void releaseHapData(HaplotypeData *data);

//reads in haplotype data and also does basic checks on integrity of format
//returns a populated HaplotypeData structure if successful
//throws an exception otherwise
HaplotypeData *readHaplotypeData(string filename);
HaplotypeData *readHaplotypeDataTPED(string filename);
HaplotypeData *readHaplotypeDataVCF(string filename);

//counts the number of "fields" in a string
//where a field is defined as a contiguous set of non whitespace
//characters and fields are delimited by whitespace
int countFields(const string &str);

#endif
