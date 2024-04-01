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

#include "selscan-pbar.h"

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
    bool flipped;
    vector <unsigned int> xors;
    vector <unsigned int> positions; // 0-based indexing of derived alleles ("1") (if flipped false)
    vector <unsigned int> positions2; // 0-based indexing of allele "2"
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
    ofstream* flog;

    //allocates the 2-d array and populated it with -9
    void initHapData(unsigned int nhaps, unsigned int nloci);
    void releaseHapData();
    /**
     * reads in haplotype data and also does basic checks on integrity of format
     * returns a populated HaplotypeData structure if successful
     * throws an exception otherwise
     */ 
    void readHapData(string filename);
    void readHapDataTPED(string filename);
    void readHapDataVCF(string filename, bool SKIP, double MAF);

    //double calcFreq(int locus, bool unphased)
    double calcFreq(int locus)
    {
        double total = 0;
        double freq = 0;

        //assuming no missing data in hapmap structure
        
        if(this->unphased){
            freq = hapEntries[locus].positions.size() + hapEntries[locus].positions2.size()*2;
            total = this->nhaps*2;
        }else{
            freq = hapEntries[locus].positions.size();
            total = this->nhaps;
        }
        

        // for (int hap = 0; hap < hapData.nhaps; hap++)
        // {
        //     if (hapData->data[hap][locus] != MISSING_CHAR)
        //     {
        //         if (hapData->data[hap][locus] == '1') freq += 1;
        //         else if (hapData->data[hap][locus] == '2') freq += 2;
                
        //         if (unphased) total+=2;
        //         else total++;
        //     }
        // }
        return (freq / total);
    }
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

    /***
     * @param query: Locus name
     * @returns locus ( in range [0 to nloci) )
    */
    int queryFound(string query)
    {
        if (locus_query_map.count(query)>0){
            return locus_query_map[query];
        }
        return -1;

        // for (int locus = 0; locus < this->nloci; locus++)
        // {
        //     if (this->locusName[locus].compare(query) == 0) return locus;
            
        //     This is only set when compiling windows binaries
        //     And I needed to use itoa for some reason under windows
        //     #ifdef PTW32_STATIC_LIB
        //     char buffer[100];
        //     itoa(mapData->physicalPos[locus], buffer, 10);
        //     if (string(buffer).compare(query) == 0) return locus;
        //     #else
        //     if (to_string(mapData->physicalPos[locus]).compare(query) == 0) return locus;
        //     #endif
        // }
        // return -1;
    }
};



class HapMap{
public:
    //int nloci;
    //int nhaps;

    HapData hapData;
    HapData hapData2;
    MapData mapData;

    bool SKIP;
    bool MAF;

    ofstream *flog;
    ofstream *fout;
    Bar *bar;

    bool loadHapMapData(param_main &params, int argc, char *argv[], ofstream *flog, ofstream *fout);
    //double calcFreq(string query);


    /*
    void doSkip(){
        if (SKIP) //prefilter all sites < MAF
        {
            cerr << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
            (*flog)  << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
            

            cerr << "Removed " << mapData->nloci - count << " low frequency variants.\n";
            (*flog) << "Removed " << mapData->nloci - count << " low frequency variants.\n";

            // delete [] freq;
            // freq = newfreq;
            // newfreq = NULL;

            // releaseHapData(hapData);
            // hapData = newHapData;
            // newHapData = NULL;

            // releaseMapData(mapData);
            // mapData = newMapData;
            // newMapData = NULL;
        }
    }
    */

    /***
     * @param query: Locus name
     * @returns locus ( in range [0 to nloci) ) after integrity check
    */
    int queryFound(string query){
        int queryLoc  = mapData.queryFound(query);
        
        if (queryLoc < 0)
        {
            cerr << "ERROR: Could not find specific locus query, " << query << ", in data.\n";
            exit(1);
        }else{
            double queryFreq = hapData.calcFreq(queryLoc);
            if (SKIP && (queryFreq < MAF || 1 - queryFreq < MAF))
            {
                cerr << "ERROR: EHH not calculated for " << query << ". MAF < " << MAF << ".\n";
                cerr << "\tIf you wish to calculate EHH at this locus either change --maf or set --keep-low-freq.\n";
                exit(1);
            }
            else if (!SKIP && (queryFreq == 0 || queryFreq == 1)){
                cerr << "ERROR: EHH not calculated for " << query << ". Frequency = " << queryFreq << ".\n";
                exit(1);
            }
            else
            {
                cerr << "Found " << query << " in data.\n";
            }
        }
        return  queryLoc;
    }
};


#endif
