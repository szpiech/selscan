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
// # include <unordered_map>
#include "selscan-pbar.h"

#include<queue>
#include <cmath>
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

class MyBitset{
    public:
    uint64_t* bits;
    int nbits;
    int nwords;
    int WORDSZ = 64;
    int num_1s = 0;

    MyBitset(int nbits){
        this->nbits = nbits;
        this->nwords = (nbits/WORDSZ) + 1; //idea: do ceil
        bits = new uint64_t[nwords];
        for (int i = 0; i < nwords; i++)
        {
            bits[i] = 0;
        }
    }

    int count_1s(){
        int sum = 0;
        for (int k = 0; k < nwords; k++) {
            sum += __builtin_popcountll ((uint64_t)bits[k]);
        }
        return sum;
    }

    vector<int> get_bits(){
        uint64_t bitset;
        vector<int> bitvec;
        for (int k = 0; k < nwords; k++) {
            bitset = bits[k];
            std::cout<<"B"<<bitset<<std::endl;
            while (bitset != 0) {
            uint64_t t = bitset & -bitset;
            int r = __builtin_ctzl(bitset);
            std::cout<<(k * WORDSZ + r) << std::endl;
            bitvec.push_back(k * WORDSZ + r); //idea: reserve 1 counts
            //callback(k * WORDSZ + r);
            bitset ^= t;
            }
        }
        // for (int i = 0; i < nbits; i++)
        // {
        //     bitvec.push_back(get_bit(i));
        // }
        return bitvec;
    }

    void set_bit(int bit){
        bits[bit/WORDSZ] |= (uint64_t) 1 << (bit % WORDSZ);
    }

    void clear_bit(int bit){
        bits[bit/WORDSZ] &= ~((uint64_t) 1 << (bit % WORDSZ));
    }

    bool get_bit(int bit){
        return (bits[bit/WORDSZ] & ((uint64_t) 1 << (bit % WORDSZ))) != 0;
    }

    void print(){
        for (int i = 0; i < nbits; i++)
        {
            cout << get_bit(i);
        }
        cout << endl;
    }

    void print_pos(){
        for (int i = 0; i < nbits; i++)
        {
            if(get_bit(i)==1)
                cout << i << " ";
        }
        cout << endl;
    }

    MyBitset operator^(const MyBitset& b) {
        MyBitset xor_bitset(this->nbits);
        for (int k = 0; k < this->nwords; k++) {
            xor_bitset.bits[k] = this->bits[k] ^ b.bits[k];
        }
        return xor_bitset;
    }

    ~MyBitset(){
        delete [] bits;
    }

};

struct HapEntry
{
    bool flipped = false;
    vector <unsigned int> xors;

    vector <unsigned int> xors1; //unphased
    vector <unsigned int> xors2; //unphased

    vector <unsigned int> positions; // 0-based indexing of derived alleles ("1") (if flipped false)
    int count1 = 0;
    int count2 = 0;
    vector <unsigned int> positions2; // 0-based indexing of allele "2"

    MyBitset* hapbitset;
    MyBitset* xorbitset;

};

struct MapEntry
{
    unsigned int physicalPos;
    double geneticPos;
    string locusName;
    string chr;
    int locId;
};


class HapData
{
public:
    struct HapEntry* hapEntries = NULL; //vector of haplotype entries
    unsigned int nloci;
    unsigned int nhaps;


    //string benchmark_flag = "XOR";
    string benchmark_flag = "BITSET";

    bool unphased;
    double MAF;
    bool SKIP;
    
    queue<int> skipQueue;

    
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
    void readHapDataVCF(string filename);


    void initParams(bool UNPHASED, bool SKIP, double MAF){
        this->unphased = UNPHASED;
        this->SKIP = SKIP;
        this->MAF = MAF;
    }


    void print(){
        for (int i = 0; i < 3; i++)
        {
            cout << "Locus: " << i << endl;
            cout << "Flipped: " << hapEntries[i].flipped << endl;
            cout << "Xors: ";
            for (int j = 0; j < hapEntries[i].xors.size(); j++)
            {
                cout << hapEntries[i].xors[j] << " ";
            }
            cout << endl;
            cout << "Positions: ";
            for (int j = 0; j < hapEntries[i].positions.size(); j++)
            {
                cout << hapEntries[i].positions[j] << " ";
            }
            cout << endl;

            // if(hapEntries[i].positions2.size()>0){
            //     cout << "Positions2: ";
            //     for (int j = 0; j < hapEntries[i].positions2.size(); j++)
            //     {
            //         cout << hapEntries[i].positions2[j] << " ";
            //     }
            //     cout << endl;
            // }
            
        }
    }
    //double calcFreq(int locus, bool unphased)
    double calcFreq(int locus)
    {
        double total = 0;
        double freq = 0;

        //assuming no missing data in hapmap structure
        
        if(this->unphased){
            freq = hapEntries[locus].count1 + hapEntries[locus].count2*2;
            //freq = hapEntries[locus].positions.size() + hapEntries[locus].positions2.size()*2;
            total = this->nhaps*2;
        }else{
            freq = hapEntries[locus].positions.size();
            if(benchmark_flag=="BITSET")
                //freq = hapEntries[locus].hapbitset->count_1s();
                freq = hapEntries[locus].hapbitset->num_1s;

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
    void readMapDataVCF(string filename, int expected_loci, queue<int>& skipQueue); //Physical positions only


    void print(){
        for (int i = 0; i < nloci; i++)
        {
            cout << "Locus: " << i << " Physical Pos: " << mapEntries[i].physicalPos << endl;
            // cout << "Genetic Pos: " << mapEntries[i].geneticPos << endl;
            // cout << "Locus Name: " << mapEntries[i].locusName << endl;
            // cout << "Chromosome: " << mapEntries[i].chr << endl;
        }
    }

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
