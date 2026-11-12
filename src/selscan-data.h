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

#include "gzstream.h"
#include "param_t.h"

#include <chrono>
#include <queue>
#include <deque>
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

// M9: one pass over a VCF file. Each record is checked and filtered, and
// its alleles are written to bitsets in the same pass; the bitsets of kept
// records go to HapData without a second read of the file. Kept loci are
// held as parallel vectors (8-24 bytes each) to keep the extra memory
// during loading small.
struct VCFRecords {
    int current_ngts = -1;
    size_t nloci_before_filtering = 0;
    size_t n_nonsnp = 0;                     // removed by --filter-vcf
    size_t nSkipDup = 0, nSkipFilter = 0, nSkipMaf = 0;
    // kept loci, in file order
    vector<MyBitset*> hap;                   // 1s (phased) or heterozygotes (unphased)
    vector<MyBitset*> xr;                    // unphased only: homozygous derived
    vector<int> derived;                     // XP only: derived-allele count
    queue<size_t> skipQueue;                 // record indices of removed records, in order
    deque<MapEntry> map;                     // only when the map is taken from this file
    // XP only, for every record: position, FNV-1a hash of "REF\tALT", removed by a per-file rule
    vector<size_t> physpos;
    vector<uint64_t> refalt_hash;
    vector<char> skipped;
    void freeAll() { for (auto b : hap) delete b; for (auto b : xr) delete b; hap.clear(); xr.clear(); }
};
class HapMap{
public:
    /// @brief this parameter is used to store the command line arguments, but it has two copies, this one is hapmap specific, the other one is "run" specific. This is required for multi-stat calculation.
    
    param_main p;
    vector<param_main> ps;
    set<string> chr_set;
    double MIN_MAF_CUTOFF = -1;

    // Loci (indices before any filtering) that every reader must skip
    // regardless of MAF: sites at a duplicated position. Built by
    // loadHapMapData before the haplotype readers run; empty for VCF input,
    // whose readers detect duplicates themselves.
    vector<char> preSkip;
    bool preSkipAt(size_t locus) const { return locus < preSkip.size() && preSkip[locus]; }
    void markDuplicatePositions(const vector<long>& positions, const string& source);

    std::unique_ptr<MapData> mapData;
    std::unique_ptr<HapData> hapData;
    std::unique_ptr<HapData> hapData2;

    ofstream* flog;
    // ofstream* fout;

    std::atomic<int> currentProcessed = 0;



    // P1: map entries of kept query-VCF loci, recorded during the haplotype pass
    deque<MapEntry> vcfMap;   // M9: deque, so growing it never copies and setFromEntries frees it as it goes
    VCFRecords readVCFRecords(const string &filename, bool xp, bool captureMap);
    void installVCFRecords(VCFRecords &rs, HapData &hapData, bool captureMap);
    void readHapDataVCFXP(string filename, string filename2, HapData& hapData, HapData& hapData2);





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
    void readHapDataTHAP(string filename, HapData& hap);
    void checkHapRow(const string &haps, const string &filename, int lineNum, bool firstRow);
    void readHapDataTPED(string filename, HapData& hap);
    void readHapDataVCF(string filename, HapData& hap);
    void readHapDataMSHap(string filename, HapData& hap);

    bool shouldSkipLocus(int number_of_1s, int number_of_2s, int nalleles); // int missing_count = 0
    void loadHapMapData();

};

#endif
