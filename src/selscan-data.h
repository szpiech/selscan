#ifndef __XP_IHH_DATA_H__
#define __XP_IHH_DATA_H__
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

const int MISSING = -9999;

struct HaplotypeData
{
  short **data;
  int nhaps;
  int nloci;
};

struct MapData
{
  int* physicalPos;
  double* geneticPos;
  string* locusName;
  int nloci;
  string chr;
};

//allocates the arrays and populates them with -9 or "--" depending on type
MapData* initMapData(int nloci);
void releaseMapData(MapData *data);

//reads in map data and also does basic checks on integrity of format
//returns a populated MapData structure if successful
//throws an exception otherwise
MapData* readMapData(string filename,int expected_loci);

//allocates the 2-d array and populated it with -9
HaplotypeData* initHaplotypeData(unsigned int nhaps, unsigned int nloci);
void releaseHapData(HaplotypeData *data);

//reads in haplotype data and also does basic checks on integrity of format
//returns a populated HaplotypeData structure if successful
//throws an exception otherwise
HaplotypeData* readHaplotypeData(string filename);

//counts the number of "fields" in a string
//where a field is defined as a contiguous set of non whitespace 
//characters and fields are delimited by whitespace
int countFields(const string &str);

#endif
