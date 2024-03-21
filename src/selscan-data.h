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
#include "param_t.h"

using namespace std;

const int MISSING = -9999;
const char MISSING_CHAR = '9';

/**
 * counts the number of "fields" in a string
 * where a field is defined as a contiguous set of non whitespace
 * characters and fields are delimited by whitespace
 *  @param str the string to count fields in
*/
int countFields(const string &str);


struct HapEntry
{
    //char **data;
    //data -> allpositions
    bool flipped;
    vector <unsigned int> xors;
    vector <unsigned int> positions; // 0-based indexing of derived alleles ("1") (if flipped false)
    vector <unsigned int> positions2; // 0-based indexing of allele "2"

    // int nhaps;
    // int nloci;
};

struct MapEntry
{
    int physicalPos;
    double geneticPos;
    string locusName;
    string chr;
};


class HapData
{
public:
    struct HapEntry* hapEntries = NULL; //vector of haplotype entries
    unsigned int nloci;
    unsigned int nhaps;
    bool unphased;

    //allocates the 2-d array and populated it with -9
    void initHapData(unsigned int nhaps, unsigned int nloci);
    void releaseHapData();
    /**
     * reads in haplotype data and also does basic checks on integrity of format
     * returns a populated HaplotypeData structure if successful
     * throws an exception otherwise
     */ 
    void readHapData(string filename, bool unphased);
    void readHapDataTPED(string filename, bool unphased);
    void readHapDataVCF(string filename, bool unphased);
};

class MapData
{
public:
    struct MapEntry* mapEntries = NULL; //vector of map entries
    unsigned int nloci;

    map<string, int> locus_query_map;

    //allocates the arrays and populates them with -9 or "--" depending on type
    void initMapData(int nloci);
    void releaseMapData();

    
    //reads in map data and also does basic checks on integrity of format
    //returns a populated MapData structure if successful
    //throws an exception otherwise
    void readMapData(string filename, int expected_loci, bool USE_PMAP);
    void readMapDataTPED(string filename, int expected_loci, int expected_haps, bool USE_PMAP);
    void readMapDataVCF(string filename, int expected_loci); //Physical positions only

};


class HapMap{
public:
    //int nloci;
    //int nhaps;
    

    HapData hapData;
    HapData hapData2;
    MapData mapData;

    bool loadHapMapData(param_t &params, int argc, char *argv[]);
    double calcFreq(int locus);
    double calcFreq(string query);


};


#endif
