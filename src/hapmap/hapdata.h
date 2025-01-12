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
#include <fstream>

using namespace std;

class MissingInfo{
    public:
        // int i = -1;
        // int j = -1;
        int length_from_core = -1;
        int num_samples = -1;
        double p0 = -1;
        double p1 = -1;
        double p2 = -1;
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

        MissingInfo(int length_from_core, int num_samples, double p0, double p1, double p2){
            // this->i = i;
            // this->j = j;
            this->length_from_core = length_from_core;
            this->num_samples = num_samples;
            this->p0 = p0;
            this->p1 = p1;
            this->p2 = p2;
            this->verdict =  '2';
            if (p0 >= p1 && p0 >= p2) this->verdict =  '0';
            if (p1 >= p0 && p1 >= p2) this->verdict =  '1';
            if( (abs(p0-p1)<0.0001 && p0>=p2) || abs(p1-p2)<0.0001  && p1>=p0 || abs(p0-p2)<0.0001 && p0>=p1){
                this->verdict = 'T';
            }
            //this->verdict = (abs(p0-p1)<0.0001)? 'T' :  this->verdict;
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
    bool MULTI_CHR = false;
    bool MULTI_MAF = false;

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

        /// @brief 
    /// @param i : hapIndex
    /// @param j : locusIndex
    /// @param length_from_core 
    /// @param num_samples 
    /// @param p0 
    /// @param p1 
    void insert_into_missing_matrix(int i, int j, int length_from_core, int num_samples, double p0, double p1, double p2){
        std::lock_guard<std::mutex> lock(map_mutex);
        {
            if(missingMatrix.find(make_pair(i,j)) == missingMatrix.end()){
                missingMatrix[make_pair(i,j)] = vector<MissingInfo>();
                
            }
            missingMatrix[make_pair(i,j)].push_back(MissingInfo(length_from_core, num_samples, p0, p1, p2));
        }
        
    }
    

    void assignVerdict(){
        lock_guard<mutex> lock(map_mutex);
        {
            if(!unphased){
                for(auto && [key, value]: missingMatrix){
                    double count1 = 0;
                    double count0 = 0;
                    for(auto && info: value){
                        //cout<<key.first<<", "<<key.second <<" "<<info.length_from_core<<" "<< info.num_samples << " "<<info.p0<<" "<<info.p1<<" "<<info.verdict<<endl;
                        //info.print(std::ref(map_mutex));
                        // if(info.verdict == '0'){
                        //     count0++;
                        // }else if(info.verdict == '1'){
                        //     count1+=;
                        // }
                        count1+=info.p1;
                        count0+=info.p0;
                        // if(count0>count1){
                        //     hapEntries[key.second].hapbitset->set_bit(key.first);
                        // }
                    }
                    if(count1>count0){
                        hapEntries[key.second].hapbitset->set_bit(key.first);
                        hapEntries[key.second].hapbitset->num_1s++;
                    }
                }
            }else{
                for(auto && [key, value]: missingMatrix){
                    double count1 = 0;
                    double count0 = 0;
                    double count2 = 0;
                    for(auto && info: value){
                        cout<<key.first<<", "<<key.second <<" "<<info.length_from_core<<" "<< info.num_samples << " p0"<<info.p0<<" p1"<<info.p1<<" p2"<<info.p2<<" v"<<info.verdict<<endl;
                        //info.print(std::ref(map_mutex));
                        // if(info.verdict == '0'){
                        //     count0++;
                        // }else if(info.verdict == '1'){
                        //     count1++;
                        // }else if(info.verdict == '2'){
                        //     count2++;
                        // }
                         count1+=info.p1;
                        count0+=info.p0;
                        count2+=info.p2;

                        // if(count1>=count0 && count1>=count2){
                        //     hapEntries[key.second].hapbitset->set_bit(key.first);
                        //     hapEntries[key.second].hapbitset->num_1s++;
                        //     //addAllele1(key.second, key.first);  
                        // }else if(count2>=count0 && count2>=count1){
                        //     //addAllele2(key.second, key.first);  
                        //     hapEntries[key.second].xorbitset->set_bit(key.first);
                        //     hapEntries[key.second].xorbitset->num_1s++;
                        // }
                    }
                    if(count1>=count0 && count1>=count2){
                            hapEntries[key.second].hapbitset->set_bit(key.first);
                            hapEntries[key.second].hapbitset->num_1s++;
                            //addAllele1(key.second, key.first);  
                        }else if(count2>=count0 && count2>=count1){
                            //addAllele2(key.second, key.first);  
                            hapEntries[key.second].xorbitset->set_bit(key.first);
                            hapEntries[key.second].xorbitset->num_1s++;
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

        //open filestream
        ofstream happrintFile("happrint.txt");
        if(!happrintFile.is_open()){
            cerr<<"ERROR: Debug log file (happrint.txt) could not be opened."<<endl;
            exit(1);
        }

        if(nloci==0){
            nloci = this->nloci;
        }

        for(int locus_after_filter = 0; locus_after_filter < nloci; locus_after_filter++){
            // cout<<locus_after_filter<<"::: ";
            // hapEntries[locus_after_filter].hapbitset->print_pos();
            happrintFile<<locus_after_filter<<"::: ";
            hapEntries[locus_after_filter].hapbitset->print_pos_to_file(&happrintFile);
        }

        happrintFile.close();

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
        //return nhaps - hapEntries[locus].xorbitset->num_1s; // why?
        return hapEntries[locus].xorbitset->num_1s;
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

    inline double get_maf(int locus){
        if(MISSING_ALLOWED){
                throw "not implemented";
            }
        if(unphased){
            throw "not implemented";
            int n_c0 = get_n_c0(locus);
            int n_c1 = get_n_c1(locus);
            int n_c2 = get_n_c2(locus);
            int n_c_missing = get_n_c_missing(locus);
            int total = n_c0 + n_c1 + n_c2 + n_c_missing;
            if(total == 0){
                return 0;
            }
            double maf = (n_c1 + 2*n_c2) / (double) (nhaps*2);
            return (maf > 0.5)? 1-maf : maf;
        }else{
            
            int n_c0 = get_n_c0(locus);
            int n_c1 = get_n_c1(locus);
            int total = n_c0 + n_c1;
            if(total == 0){
                return 0;
            }
            double maf = (double) n_c1 / nhaps;
            
            if (maf > 0.5){
                //cout<<n_c0<< " "<<n_c1<<" "<< 1-maf<<endl;
                return 1-maf;
            } else{
                //cout<<n_c0<< " "<<n_c1<<" "<< maf<<endl;
                return maf;
            }
        }
        
    }

    inline int get_n_c_missing(int locus)
    {
        if(unphased){
            return hapEntries[locus].missbitset->num_1s;
        }
        return hapEntries[locus].xorbitset->num_1s;
    }

    void addAlleleMissing(int locus, int hapId){
        if(unphased){
            hapEntries[locus].missbitset->set_bit(hapId);
            hapEntries[locus].missbitset->num_1s+= 1;
        }else{
            hapEntries[locus].xorbitset->set_bit(hapId);
            hapEntries[locus].xorbitset->num_1s+= 1;
        }
        
    }

    void addAllele1(int locus, int hapId){
        hapEntries[locus].hapbitset->set_bit(hapId);
        hapEntries[locus].hapbitset->num_1s+= 1;
    }

    void addAllele2(int locus, int hapId){
        if(!unphased){
            throw "ERROR: not unphased";
            exit(1);
        }
        hapEntries[locus].xorbitset->set_bit(hapId);
        hapEntries[locus].xorbitset->num_1s+= 1;
    }

    private:
        bool INIT_SUCCESS = false;

};


#endif