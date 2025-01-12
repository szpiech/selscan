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

#include "selscan-cli.h"

#include <sstream>
#include <algorithm>

#include <string>
#include <iostream>
#include <fstream>
#include <omp.h>
#include "gzstream.h"
#include "param_t.h"

#include <chrono>
#include <queue>
#include <cmath>

#include "hapmap/bitset.h"
#include "hapmap/hapdata.h"
#include "hapmap/mapdata.h"

#include <memory>
#include <atomic>
#include <set>

// #include "filetype/vcf.h"
// #include "filetype/vcf_serial.h"
// #include "filetype/hap.h"
// #include "filetype/hap_serial.h"

using namespace std;

const int MISSING = -9999;  // TO_BE_DELETED
const char MISSING_CHAR = '9';  // TO_BE_DELETED

class HapMap{
public:
    /// @brief this parameter is used to store the command line arguments, but it has two copies, this one is hapmap specific, the other one is "run" specific. This is required for multi-stat calculation.
    
    param_main p;
    vector<param_main> ps;
    set<string> chr_set;
    double MIN_MAF_CUTOFF = -1;

    std::unique_ptr<MapData> mapData;
    std::unique_ptr<HapData> hapData;
    std::unique_ptr<HapData> hapData2;

    ofstream* flog;
    // ofstream* fout;

    std::atomic<int> currentProcessed = 0;

    char randomZeroOrOne() {
        // Create a random device to seed the generator
        std::random_device rd; 
        
        // Use a Mersenne Twister engine
        std::mt19937 gen(rd()); 
        
        // Define a uniform distribution for integers between 0 and 1
        std::uniform_int_distribution<int> dist(0, 1);
        
        // Return the random number
        int randomNum  = dist(gen);

         // Convert integer 0 or 1 to character '0' or '1'
        char randomChar = static_cast<char>(randomNum + '0');
        return randomChar;
    }

    pair<int, int> countFieldsAndOnes(const string &str)
    {
        string::const_iterator it;
        int ones = 0;
        int result;
        int numFields = 0;
        int seenChar = 0;
        for (it = str.begin() ; it < str.end(); it++)
        {
            if(*it == '1'){
                ones++;
            }
            result = isspace(*it);
            if (result == 0 && seenChar == 0)
            {
                numFields++;
                seenChar = 1;
            }
            else if (result != 0)
            {
                seenChar = 0;
            }
        }
        return make_pair(numFields, ones);
    }

    int countFields(const string &str)
    {
        string::const_iterator it;
        int result;
        int numFields = 0;
        int seenChar = 0;
        for (it = str.begin() ; it < str.end(); it++)
        {
            result = isspace(*it);
            if (result == 0 && seenChar == 0)
            {
                numFields++;
                seenChar = 1;
            }
            else if (result != 0)
            {
                seenChar = 0;
            }
        }
        return numFields;
    }

    HapMap(param_main &p, ofstream* flog){
        this->p = p;
        this->flog = flog;
        // this->fout = fout;

    }

    HapMap(param_main &p, vector<param_main>&ps, ofstream* flog){
        this->p = p;
        this->flog = flog;
        this->ps = ps;
        // this->fout = fout;

    }

    ~HapMap(){
    }

    /**
     * reads in haplotype data and also does basic checks on integrity of format
     * returns a populated HaplotypeData structure if successful
     * throws an exception otherwise
     */ 
    void initParamsInHap(HapData& hap);
    void readHapData(string filename, HapData& hap);
    void readHapDataTPED(string filename, HapData& hap);
    void readHapDataVCF(string filename, HapData& hap);
    void readHapDataVCFMissing(string filename, HapData& hap);
    void readHapDataMSHap(string filename, HapData& hap);

    bool skipThisLocus(int number_of_1s, int number_of_2s, int nalleles, int missing_count = 0);
    void loadHapMapData();

};

#endif
