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

#ifndef __DATA_H__
#define __DATA_H__


#include <string>
#include <iostream>
#include <fstream>
#include <omp.h>
#include "gzstream.h"
#include "param_t.h"
// # include <unordered_map>
#include "selscan-pbar.h"

#include<queue>
#include <cmath>

#include "hapmap/bitset.h"
#include "hapmap/hapdata.h"
#include "hapmap/mapdata.h"
//#include "hapmap/hapmap.h"


#include "filetype/data_reader.h"
#include "filetype/vcf.h"


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
pair<int, int> countFieldsAndOnes(const string &str);



// int countFields(const string &str)
// {
//     string::const_iterator it;
//     int result;
//     int numFields = 0;
//     int seenChar = 0;
//     for (it = str.begin() ; it < str.end(); it++)
//     {
//         result = isspace(*it);
//         if (result == 0 && seenChar == 0)
//         {
//             numFields++;
//             seenChar = 1;
//         }
//         else if (result != 0)
//         {
//             seenChar = 0;
//         }
//     }
//     return numFields;
// }

// pair<int, int> countFieldsAndOnes(const string &str)
// {
//     string::const_iterator it;
//     int ones = 0;
//     int result;
//     int numFields = 0;
//     int seenChar = 0;
//     for (it = str.begin() ; it < str.end(); it++)
//     {
//         if(*it == '1'){
//             ones++;
//         }
//         result = isspace(*it);
//         if (result == 0 && seenChar == 0)
//         {
//             numFields++;
//             seenChar = 1;
//         }
//         else if (result != 0)
//         {
//             seenChar = 0;
//         }
//     }
//     return make_pair(numFields, ones);
// }


class HapMap{
public:
    //int nloci;
    //int nhaps;

    HapData hapData;
    HapData hapData2;
    MapData mapData;

    ofstream *flog;
    ofstream *fout;
    Bar *bar;

    bool loadHapMapData(param_main &params, int argc, char *argv[], ofstream *flog, ofstream *fout);
    //double calcFreq(string query);
    bool is_gzipped(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        
        if (!file.is_open()) {
            std::cerr << "Error: Unable to open file " << filename << std::endl;
            return false;
        }

        // Read the first two bytes
        std::vector<unsigned char> buffer(2);
        file.read(reinterpret_cast<char*>(buffer.data()), buffer.size());

        // Close the file
        file.close();

        // Check if the first two bytes are the gzip magic numbers
        return buffer[0] == 0x1F && buffer[1] == 0x8B;
    }
    void handleData(string filename){
        DataReader* dataReader;
        string flag="--vcf";

        if(flag=="--vcf"){
            if(is_gzipped(filename)){
                dataReader = new VCFSerialReader(filename, hapData);
            }else{
                dataReader = new VCFParallelReader(filename, hapData);
            }
        }
        
        dataReader->get_nloci_nhaps(); 
        int* num1s_per_loci = new int[hapData.nloci];
        int* num2s_per_loci = new int[hapData.nloci];
        queue<int> skiplist;
        //dataReader->n1s_n2s_q(num1s_per_loci, num2s_per_loci, skiplist);
        delete[] num1s_per_loci;
        delete[] num2s_per_loci;
    }

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
             cerr << "Found " << query << " in data.\n";
            // double queryFreq = hapData.calcFreq(queryLoc);
            // if (hapData.SKIP && (queryFreq < hapData.MAF || 1 - queryFreq < hapData.MAF))
            // {
            //     cerr << "ERROR: EHH not calculated for " << query << ". MAF < " << hapData.MAF << ".\n";
            //     cerr << "\tIf you wish to calculate EHH at this locus either change --maf or set --keep-low-freq.\n";
            //     exit(1);
            // }
            // else if (!hapData.SKIP && (queryFreq == 0 || queryFreq == 1)){
            //     cerr << "ERROR: EHH not calculated for " << query << ". Frequency = " << queryFreq << ".\n";
            //     exit(1);
            // }
            // else
            // {
            //     cerr << "Found " << query << " in data.\n";
            // }
        }
        return  queryLoc;
    }
};


#endif
