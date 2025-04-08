#ifndef __IHH12_H__
#define __IHH12_H__

#include "selscan-stats.h"
#include "../thread_pool.h"

#include <thread>
#include <unordered_map>

#define ACTION_ON_ALL_SET_BITS(hapbitset, ACTION)         \
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

using namespace std;

class IHH12_ehh_data{
    public:
    double prev_ehh12_before_norm = twice_num_pair(nhaps);
    double curr_ehh12_before_norm = 0;

    // double prev_ehh12d1_before_norm = 0;
    // double curr_ehh12d1_before_norm = 0;

    double prev_ehh_before_norm = twice_num_pair(nhaps);
    double curr_ehh_before_norm = 0;

    int n_c[2] = {0,0};
    int nhaps;
    int totgc = 0;
    int* group_count;
    bool* isDerived;
    int* group_id;

    MyBitset* v;

    ~IHH12_ehh_data(){
        delete[] group_count;
        delete[] isDerived;
        delete[] group_id;
    }

    inline double square_alt(int n){
        return n*n;
    }

    inline double twice_num_pair(int n){
        if(n < 2){
            return 0;
        }
        // return 2*nCk(n, 2);
        return 2*nCk(n, 2);
    }
    void init(int nhaps, MyBitset* positions ){
        this->nhaps = nhaps;
        group_count = new int[nhaps];
        isDerived = new bool[nhaps];
        group_id = new int[nhaps];

        //will be vectorized with compile time flags
        for(int i = 0; i< nhaps; i++){
            group_count[i] = 0;
            group_id[i] = 0;
            isDerived[i] = false;
        }
        v = positions;

        //updated in pooled
        n_c[0] = nhaps - v->num_1s;
        n_c[1] =  v->num_1s;
    }

    void initialize_core(bool ALT){
        if(n_c[1] == 0 || n_c[0] == 0){
            cout<<"WARNING: cannot use this as core"<<endl;
        }

        group_count[1] = n_c[1];
        group_count[0] = n_c[0];
        totgc+=2;

        ACTION_ON_ALL_SET_BITS(v, {
            isDerived[set_bit_pos] = true;
            group_id[set_bit_pos] = 1;
        });

        curr_ehh_before_norm = (ALT?   square_alt(n_c[0]) +  square_alt(n_c[1]) :  twice_num_pair(n_c[0]) + twice_num_pair(n_c[1])); 

        double firstFreq_before_norm = (n_c[1] > 1) ? twice_num_pair(n_c[1]) : 0;
        double secondFreq_before_norm =(n_c[0] > 1) ?  twice_num_pair(n_c[0]): 0;
        double comboFreq_before_norm = (nhaps > 1) ? twice_num_pair(nhaps) : 0;
        
        //double normfac = twice_num_pair(nhaps);
        curr_ehh12_before_norm = curr_ehh_before_norm  - firstFreq_before_norm - secondFreq_before_norm + comboFreq_before_norm;
       //curr_ehh12_before_norm = curr_ehh_before_norm + 2.0*top1*top2; 
        //DEBUG:: cout<<"t1 t2 "<<n_c[1]<<" "<<n_c[0]<<" "<<firstFreq/normfac<<" "<<secondFreq/normfac<<" "<< comboFreq/normfac<<endl; 
    }

}; 

class IHH12 : public SelscanStats{
    public:
        IHH12(const std::unique_ptr<HapMap>&  hm, param_main& params) : SelscanStats(hm, params){
            max_extend = ( p.MAX_EXTEND <= 0) ? hm->mapData->mapEntries[hm->mapData->nloci-1].physicalPos -  hm->mapData->mapEntries[0].physicalPos : p.MAX_EXTEND  ;
        }
        void main();
        void findMaxTwo(int* arr, int n, int &max1, int &max2);
        //void findMaxK(int* arr, int n, int &max1, int &max2, int k);
        std::mutex mutex_log; //static pthread_mutex_t mutex_log;

    protected:
        int max_extend;
        //void static thread_main(int tid, unordered_map<int, vector<int> >& m, IHH12* obj);
        void updateEHH_from_split(const unordered_map<int, vector<int> > & m, IHH12_ehh_data* p);
    
    private:
        double calc_ihh12_at_locus(int locus);
        double calc_ehh_unidirection(int locus, bool downstream);
};





#endif