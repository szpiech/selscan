#ifndef __HAPDATA_H__
#define __HAPDATA_H__

#include <vector>
#include <string>
#include "bitset.h"
#include <queue>
#include <iostream>




using namespace std;

struct HapEntry
{
    bool flipped = false;
    vector <unsigned int> xors;

    vector <unsigned int> xors1; //unphased
    vector <unsigned int> xors2; //unphased

    vector<unsigned int> g[3];

    vector <unsigned int> positions; // 0-based indexing of derived alleles ("1") (if flipped false)
    int count1 = 0;
    int count2 = 0;
    vector <unsigned int> positions2; // 0-based indexing of allele "2"

    MyBitset* hapbitset;
    MyBitset* xorbitset;
};



class HapData
{
public:
    void do_xor(); //experimental

    struct HapEntry* hapEntries = NULL; //vector of haplotype entries
    int nloci;
    int nhaps;

    string benchmark_flag = "XOR";
    string benchmark_flag2 = ""; //"FLIP";

    //string benchmark_flag = "BITSET";
    //string benchmark_flag = "FLIP_ONLY";
    //string benchmark_flag = "BASIC";

    string DEBUG_FLAG = "VCF";
    //string DEBUG_FLAG = "";

    bool unphased;
    double MAF;
    bool SKIP;
    int num_threads;
    bool LOW_MEM = false;
    
    queue<int> skipQueue;
    
    ofstream* flog;

    ~HapData(){
        if(LOW_MEM){
            releaseHapData_bitset();
        }else{
            releaseHapData();
        }
    }

    //allocates the 2-d array and populated it with -9,
    /** Sets up structure according to nhaps and nloci
     * 
    */
    void initHapData(int nhaps, int nloci);
    void releaseHapData();
    /**
     * reads in haplotype data and also does basic checks on integrity of format
     * returns a populated HaplotypeData structure if successful
     * throws an exception otherwise
     */ 
    void readHapData(string filename);
    void readHapDataTPED(string filename);
    void readHapDataVCF(string filename);

    void initHapData_bitset(int nhaps, unsigned int nloci);
    void releaseHapData_bitset();
    void readHapDataVCF_bitset(string filename);
    void readHapData_bitset(string filename);

    void initParams(bool UNPHASED, bool SKIP, double MAF, int num_threads, ofstream* flog, bool LOW_MEM){
        this->unphased = UNPHASED;
        this->SKIP = SKIP;
        this->MAF = MAF;
        this->num_threads = num_threads;
        this->flog = flog;
        this->LOW_MEM = LOW_MEM;
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


    void print(){
        for(int locus_after_filter = 0; locus_after_filter < this->nloci; locus_after_filter++){
            cout<<locus_after_filter<<"::: ";
            for(int i = 0; i < hapEntries[locus_after_filter].positions.size(); i++){
                cout<<hapEntries[locus_after_filter].positions[i]<<" ";
            }
            cout<<endl;
            
        }

        /*
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
        */
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

            if(LOW_MEM){
                //freq = hapEntries[locus].hapbitset->count_1s();
                freq = hapEntries[locus].hapbitset->num_1s + (hapEntries[locus].xorbitset->num_1s-hapEntries[locus].hapbitset->num_1s)*2;
            }
        }else{
            freq = hapEntries[locus].positions.size();
            if(LOW_MEM){
                //freq = hapEntries[locus].hapbitset->count_1s();
                freq = hapEntries[locus].hapbitset->num_1s;
            }
                

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

    inline unsigned int get_n_c0(int locus)
    {
        unsigned int non_flipped_1s = 0;
        non_flipped_1s = LOW_MEM? hapEntries[locus].hapbitset->num_1s: hapEntries[locus].positions.size();
        unsigned int result = hapEntries[locus].flipped ? non_flipped_1s : nhaps - non_flipped_1s;
        return result;
        //return hapEntries[locus].count1;
    }
    inline unsigned int get_n_c1(int locus)
    {
        unsigned int non_flipped_1s = 0;
        non_flipped_1s = LOW_MEM? hapEntries[locus].hapbitset->num_1s: hapEntries[locus].positions.size();
        unsigned int result = hapEntries[locus].flipped ? nhaps - non_flipped_1s : non_flipped_1s;
        return result;
        //return hapEntries[locus].count1;
    }


    template<typename T>
    void getThreeUnphasedGroups(const std::vector<T>& new1, const std::vector<T>& old1, const std::vector<T>& new2, const std::vector<T>& old2, std::vector<T> g[]) {
        std::vector<T> symdif1;
        std::vector<bool> symdif1_isOld;
        std::vector<T> symdif2;
        std::vector<bool> symdif2_isOld;

        auto it1 = new1.begin();
        auto it2 = old1.begin();
        while (it1 != new1.end() && it2 != old1.end()) {
            if (*it1 < *it2) {
                symdif1.push_back(*it1);
                ++it1;
                symdif1_isOld.push_back(false);
                //in1 not in2  (new1)  
            } else if (*it2 < *it1) {
                symdif1.push_back(*it2);
                ++it2;
                symdif1_isOld.push_back(true);
                //in2 not in1   (old1)
            } else {
                // If the elements are equal, skip them in both vectors
                ++it1;
                ++it2;
                //in both
            }
        }

        while (it1 != new1.end()) { // Copy any remaining elements from vec1
            symdif1.push_back(*it1);
            ++it1;
            symdif1_isOld.push_back(false);

            //in1 not in2
        }
        while (it2 != old1.end()) {  // Copy any remaining elements from vec2
            symdif1.push_back(*it2);
            ++it2;
            symdif1_isOld.push_back(true);

            //in2 not in1
        }


        it1 = new2.begin();
        it2 = old2.begin();
        while (it1 != new2.end() && it2 != old2.end()) {
            if (*it1 < *it2) {
                symdif2.push_back(*it1);
                ++it1;
                symdif2_isOld.push_back(false);
                //in1 not in2  (new1)  
            } else if (*it2 < *it1) {
                symdif2.push_back(*it2);
                ++it2;
                symdif2_isOld.push_back(true);
                //in2 not in1   (old1)
            } else {
                // If the elements are equal, skip them in both vectors
                ++it1;
                ++it2;
                //in both
            }
        }
        while (it1 != new2.end()) { // Copy any remaining elements from vec1
            symdif2.push_back(*it1);
            ++it1;
            symdif2_isOld.push_back(false);
            //in1 not in2
        }
        while (it2 != old2.end()) {  // Copy any remaining elements from vec2
            symdif2.push_back(*it2);
            ++it2;
            symdif2_isOld.push_back(true);
            //in2 not in1
        }


        it1 = symdif1.begin();
        it2 = symdif2.begin();

        auto it1_bool = symdif1_isOld.begin();
        auto it2_bool = symdif2_isOld.begin();

        //std::vector<T> g[3];

        int txor[3][3]; //old->new
        txor[0][0] = 0;
        txor[0][1] = 2;
        txor[0][2] = 1;
        txor[1][0] = 1;
        txor[1][1] = 0;
        txor[1][2] = 2;
        txor[2][0] = 2;
        txor[2][1] = 1;
        txor[2][2] = 0;

        while (it1 != symdif1.end() && it2 != symdif2.end()) {
            if (*it1 < *it2) {
                //not in 2
                if((*it1_bool)==true){ //old 1
                    //new 0
                    g[txor[1][0]].push_back(*it1);
                }else{ //new 1
                    //old 0
                    g[txor[0][1]].push_back(*it1);
                }
                ++it1;
                ++it1_bool;
                //1-2  (new1)  
            } else if (*it2 < *it1) {
                //not in 1
                if((*it2_bool)==true){ //old 2
                    //new 0
                    g[txor[2][0]].push_back(*it2);
                }else{ //new 2
                    //old 0
                    g[txor[0][2]].push_back(*it2);
                }
                ++it2;
                ++it2_bool;
                //2-1   (old1)
            } else {
                // If the elements are equal, skip them in both vectors
                if((*it1_bool)==true){ //old 1
                    //new 2
                    g[txor[1][2]].push_back(*it1);
                }else{ //new 1
                    //old 2
                    g[txor[2][1]].push_back(*it1);
                }
                ++it1;
                ++it2;
                ++it1_bool;
                ++it2_bool;
                //1 and 2
            }
        }
        while (it1 != symdif1.end()) { // Copy any remaining elements from vec1
            //not in 2
            if((*it1_bool)==true){ //old 1
                //new 0
                g[txor[1][0]].push_back(*it1);
            }else{ //new 1
                //old 0
                g[txor[0][1]].push_back(*it1);
            }
            ++it1;
            ++it1_bool;
            //1-2
        }
        while (it2 != symdif2.end()) {  // Copy any remaining elements from vec2
            //not in 1
            if((*it2_bool)==true){ //old 2
                //new 0
                g[txor[2][0]].push_back(*it2);
            }else{ //new 2
                //old 0
                g[txor[0][2]].push_back(*it2);
            }
            ++it2;
            ++it2_bool;
            //2-1
        }
    }

    private:
        bool INIT_SUCCESS = false;

};


#endif