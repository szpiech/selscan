#ifndef __HAPDATA_H__
#define __HAPDATA_H__

#include <vector>
#include <string>
#include "bitset.h"
#include <queue>
#include <iostream>
#include <random>
#include<map>

#include<mutex>
#include <thread>
#include<optional>

using namespace std;

class MissingInfo{
    public:
        // int i = -1;
        // int j = -1;
        int length_from_core = -1;
        int num_samples = -1;
        double p0 = -1;
        double p1 = -1;
        char verdict = 'N'; //0,1,2 for tie

        MissingInfo(int length_from_core, int num_samples, double p0, double p1){
            // this->i = i;
            // this->j = j;
            this->length_from_core = length_from_core;
            this->num_samples = num_samples;
            this->p0 = p0;
            this->p1 = p1;
            this->verdict = (p0 > p1)? '0' : '1';
            this->verdict = (abs(p0-p1)<0.0001)? 'T' :  this->verdict;
        }
        void print(std::mutex& map_mutex){
            std::lock_guard<std::mutex> lock(map_mutex);
            cout<<length_from_core<<" "<< num_samples << " "<<p0<<" "<<p1<<" "<<verdict<<endl;
        }

};

#define FACTION_ON_ALL_SET_BITS(hapbitset, ACTION)         \
    for (int k = 0; k < (hapbitset->nwords); k++) {             \
        uint64_t bitset = (hapbitset->bits)[k];                 \
        while (bitset != 0) {                        \
            uint64_t t = bitset & -bitset;           \
            int r = __builtin_ctzl(bitset);          \
            int set_bit_pos = (k * 64 + r);          \
            bitset ^= t;                             \
            ACTION;                                  \
        }                                            \
    }


struct HapEntry
{
    bool flipped = false;
    
    MyBitset* hapbitset; // stores 1s (set bits indicate 1s)
    MyBitset* xorbitset; // for unphased stores 2s, for phased stores consecutive loci XORs 

    //for missing
    MyBitset* missbitset; // stores missing bits
};

class HapData
{
public:
    map< pair<int, int> , vector<MissingInfo> > missingMatrix;
    std::mutex map_mutex;  

    void xor_for_phased_and_unphased(); //experimental

    struct HapEntry* hapEntries = NULL; //vector of haplotype entries
    int nloci;
    int nhaps;
    bool unphased = false;
    bool MISSING_ALLOWED = false;

    // string MULTI_CHR_LIST = "";
    // double MAF;
    // bool SKIP;
    // int num_threads;
    // bool LOW_MEM = true;

    queue<int> skipQueue;
    
    //ofstream* flog;

    ~HapData();

    /// @brief 
    /// @param i : hapIndex
    /// @param j : locusIndex
    /// @param length_from_core 
    /// @param num_samples 
    /// @param p0 
    /// @param p1 
    void insert_into_missing_matrix(int i, int j, int length_from_core, int num_samples, double p0, double p1){
        std::lock_guard<std::mutex> lock(map_mutex);
        {
            if(missingMatrix.find(make_pair(i,j)) == missingMatrix.end()){
                missingMatrix[make_pair(i,j)] = vector<MissingInfo>();
                
            }
            missingMatrix[make_pair(i,j)].push_back(MissingInfo(length_from_core, num_samples, p0, p1));
        }
        
    }

    void assignVerdict(){
        lock_guard<mutex> lock(map_mutex);
        {
            for(auto && [key, value]: missingMatrix){
                int count1 = 0;
                int count0 = 0;
                for(auto && info: value){
                    //cout<<key.first<<", "<<key.second <<" "<<info.length_from_core<<" "<< info.num_samples << " "<<info.p0<<" "<<info.p1<<" "<<info.verdict<<endl;
                    //info.print(std::ref(map_mutex));
                    if(info.verdict == '0'){
                        count0++;
                    }else if(info.verdict == '1'){
                        count1++;
                    }
                    if(count0>count1){
                        hapEntries[key.second].hapbitset->set_bit(key.first);
                    }
                }
            }
        }
    }

    void printMissingMatrix(){
        lock_guard<mutex> lock(map_mutex);
        {
            for(auto && [key, value]: missingMatrix){
                int count1 = 0;
                int count0 = 0;
                for(auto && info: value){
                    cout<<key.first<<", "<<key.second <<" "<<info.length_from_core<<" "<< info.num_samples << " "<<info.p0<<" "<<info.p1<<" "<<info.verdict<<endl;
                }
            }
        }
        
    }

    //allocates the 2-d array and populated it with -9,
    /** Sets up structure according to nhaps and nloci
     * 
    */
    void initHapData(int nhaps, int nloci);
    void releaseHapData();
    


    // void readHapDataVCF_bitset(string filename);
    // void readHapData_bitset(string filename);


    void print(int nloci=0){
        if(nloci==0){
            nloci = this->nloci;
        }

        for(int locus_after_filter = 0; locus_after_filter < nloci; locus_after_filter++){
            // cout<<locus_after_filter<<"::: ";
            // hapEntries[locus_after_filter].hapbitset->print_pos();
            cout<<locus_after_filter<<"::: ";
            hapEntries[locus_after_filter].hapbitset->print_pos();
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

    void print_by_ij(){
        cout<<"Printing by ij"<<endl;
        for(int locus_after_filter = 0; locus_after_filter < this->nloci; locus_after_filter++){
            FACTION_ON_ALL_SET_BITS(hapEntries[locus_after_filter].hapbitset, {
                 cout<< set_bit_pos << " " << locus_after_filter<< " " <<" "<<1 <<endl;
            });
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
    

    /// @brief assumes no missing
    /// @param locus 
    /// @return 
    double calcFreq(int locus)
    {
        //assuming no flip
        double total = 0;
        double freq = 0;
        if(this->unphased){
            freq = hapEntries[locus].hapbitset->num_1s + (hapEntries[locus].xorbitset->num_1s)*2;
            total = this->nhaps*2;
        }else{
            freq = hapEntries[locus].hapbitset->num_1s;
            total = this->nhaps;
        }
        return (freq / total);
    }

    double calcMissingFreq(int locus){
        if(!MISSING_ALLOWED){
            throw "missing not allowed";
        }
        return hapEntries[locus].xorbitset->num_1s / (double) nhaps;
    }

    inline int get_n_c2(int locus)
    {
        return nhaps - hapEntries[locus].xorbitset->num_1s;
        // // if flip was enabled
        // int non_flipped_1s = 0;
        // non_flipped_1s = LOW_MEM? hapEntries[locus].hapbitset->num_1s: hapEntries[locus].positions.size();
        // int result = hapEntries[locus].flipped ? non_flipped_1s : nhaps - non_flipped_1s;
        // return result;
    }


    inline int get_n_c0(int locus)
    {
        if(MISSING_ALLOWED){
            return nhaps - hapEntries[locus].hapbitset->num_1s - hapEntries[locus].xorbitset->num_1s; //hapEntries[locus].xorbitset->num_1s is missing count
            throw "not implemented";
        }else{
            return nhaps - hapEntries[locus].hapbitset->num_1s;

        }
        // // if flip was enabled
        // int non_flipped_1s = 0;
        // non_flipped_1s = LOW_MEM? hapEntries[locus].hapbitset->num_1s: hapEntries[locus].positions.size();
        // int result = hapEntries[locus].flipped ? non_flipped_1s : nhaps - non_flipped_1s;
        // return result;
    }

    inline int get_n_c1(int locus)
    {
        return  hapEntries[locus].hapbitset->num_1s;
        // // if flip was enabled
        // int non_flipped_1s = 0;
        // non_flipped_1s = LOW_MEM? hapEntries[locus].hapbitset->num_1s: hapEntries[locus].positions.size();
        // int result = hapEntries[locus].flipped ? nhaps - non_flipped_1s : non_flipped_1s;
        // return result;
    }

    inline int get_n_c_missing(int locus)
    {
        if(unphased){
            return hapEntries[locus].missbitset->num_1s;
        }
        return hapEntries[locus].xorbitset->num_1s;
    }

    private:
        bool INIT_SUCCESS = false;

};


#endif