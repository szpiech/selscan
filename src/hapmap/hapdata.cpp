#include "hapdata.h"
#include "../gzstream.h"
#include "../selscan-cli.h"
#include <algorithm>
#include <sstream>
#include <memory>
#include <limits>
#include<set>

HapData::~HapData(){
    releaseHapData();
}

///
/// @brief Sets up structure according to nhaps and nloci
void HapData::initHapData(int nhaps, int nloci)
{
    if (nhaps > std::numeric_limits<int>::max() || nloci > std::numeric_limits<int>::max())
    {
        cerr << "ERROR: number of haplotypes (" << nhaps << ") and number of loci (" << nloci << ") exceed maximum value.\n";
        throw 0;
    }

    if (nhaps < 1 || nloci < 1)
    {
        cerr << "ERROR: number of haplotypes (" << nhaps << ") and number of loci (" << nloci << ") must be positive.\n";
        throw 0;
    }
    cout << "=====" << endl;
    std::cout << "After filtering: \n"
          << "  Haplotypes: " << nhaps << "\n"
          << "  Loci:       " << nloci << std::endl;
    cout << "=====" << endl;
    this->hapEntries = new struct HapEntry[nloci];
    this->nhaps = nhaps;
    this->nloci = nloci;
    
    for (int j = 0; j < nloci; j++){
        hapEntries[j].hapbitset = new MyBitset(nhaps);
        hapEntries[j].xorbitset = new MyBitset(nhaps);
        if(MISSING_ALLOWED && unphased){
            hapEntries[j].missbitset = new MyBitset(nhaps);
        }
    }
    INIT_SUCCESS = true;
}

void HapData::releaseHapData()
{
    if (hapEntries == NULL) return;

    //we have a MyBitset for every locus
    for (int j = 0; j < nloci; j++){
        delete hapEntries[j].hapbitset ; //MyBitset destructor called
        delete hapEntries[j].xorbitset; //MyBitset destructor called
    }
    delete [] hapEntries;
    hapEntries = NULL;
    this->nhaps = -9;
    this->nloci = -9;
    return;
}

void HapData::xor_for_phased_and_unphased(){
    // disabling XOR for now
    if(!unphased){
        if(MISSING_ALLOWED){
            return;
        }
        for(int nloci_after_filtering = 0; nloci_after_filtering < this->nloci; nloci_after_filtering++){
            if(nloci_after_filtering==0){
                // if(nloci<=1){
                //     throw "ERROR: Dataset has only 1 locus, XOR out of bound.";    
                // }
                // CHANGEXOR
                MyBitset* b1 =(hapEntries[0].hapbitset);
                MyBitset* b2 = (hapEntries[1].hapbitset);
                int sum = 0;
                for (int k = 0; k < b1->nwords; k++) {
                    hapEntries[0].xorbitset->bits[k] = b1->bits[k] ^ b2->bits[k];
                    sum += __builtin_popcountll(hapEntries[0].xorbitset->bits[k]);
                }
                hapEntries[0].xorbitset->num_1s = sum;
                
            }else{
                
                MyBitset* b1 =(hapEntries[nloci_after_filtering].hapbitset);
                MyBitset* b2 = (hapEntries[nloci_after_filtering-1].hapbitset);

                int sum = 0;
                for (int k = 0; k < b1->nwords; k++) {
                    hapEntries[nloci_after_filtering].xorbitset->bits[k] = b1->bits[k] ^ b2->bits[k];
                    sum += __builtin_popcountll(hapEntries[nloci_after_filtering].xorbitset->bits[k]);
                }
                hapEntries[nloci_after_filtering].xorbitset->num_1s = sum;
                                 
            }
        }
    }    

    // FLIP TEST
    // PHASE 4:  FLIP
    // for (int locus = 0; locus < nloci_after_filtering; locus++){
    //     if(hapEntries[locus].flipped){
    //         vector<int> zero_positions(nhaps - hapEntries[locus].positions.size());
    //         int j = 0;
    //         int front_one = hapEntries[locus].positions[j++];
    //         for(int i=0; i<nhaps; i++){
    //             if(i==front_one){
    //                 front_one = hapEntries[locus].positions[j++];
    //             }else{
    //                 zero_positions.push_back(i);
    //             }   
    //         }
    //         hapEntries[locus].positions = zero_positions;
    //     }
    // }
    // //PHASE 4: FLIP

    /*
    if(benchmark_flag == "XOR" || benchmark_flag == "FLIP_ONLY" ){
        for (int locus_after_filter = 0; locus_after_filter < this->nloci; locus_after_filter++){
            if(hapEntries[locus_after_filter].positions.size() > this->nhaps/2){
                hapEntries[locus_after_filter].flipped = true;

                vector<int> copy_pos;
                copy_pos.reserve(this->nhaps - hapEntries[locus_after_filter].positions.size());
                int cnt = 0;
                for(int i = 0; i< this->nhaps; i++){
                    int curr =  hapEntries[locus_after_filter].positions[cnt];
                    if(i==curr){
                        cnt++;
                    }else{
                        copy_pos.push_back(i);
                    }
                }
                
                this->hapEntries[locus_after_filter].positions = copy_pos;
                copy_pos.clear();
            }else{
                this->hapEntries[locus_after_filter].flipped = false;
            }
        }
    }
    */
}

