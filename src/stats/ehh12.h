#include<unordered_map>

#ifndef __SELSCAN_EHH12_H__
#define __SELSCAN_EHH12_H__
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
#include "selscan-stats.h"
class EHH12_ehh_data{
    public:
    double curr_ehh12_before_norm = 0;
    double curr_ehh12d1_before_norm = 0;
    double curr_ehh_before_norm = 0;

    int n_c[2] = {0,0};
    int nhaps;
    int totgc = 0;
    int* group_count;
    bool* isDerived;
    int* group_id;

    MyBitset* v;

    ~EHH12_ehh_data(){
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
        curr_ehh12d1_before_norm = curr_ehh_before_norm - firstFreq_before_norm ;
       //curr_ehh12_before_norm = curr_ehh_before_norm + 2.0*top1*top2; 
        //DEBUG:: cout<<"t1 t2 "<<n_c[1]<<" "<<n_c[0]<<" "<<firstFreq/normfac<<" "<<secondFreq/normfac<<" "<< comboFreq/normfac<<endl; 
    }

}; 
class EHH12Output{
    public:
        bool print = false;
        double gdist;
        int pdist; 
        double ehh1 = 0;
        double ehh12 = 0;
        double ehh2d1 = 0;
};

class EHH12 : public SelscanStats{
    public:
        EHH12(const std::unique_ptr<HapMap>&  hm, param_main& params) : SelscanStats(hm, params){   }
        void main(string query);
        
    private:
        vector<EHH12Output> output;
        map<string, int> locus_query_map;
        void calc_ehh12_unidirection(int locus, bool downstream);
        void updateEHH_from_split(const unordered_map<int, vector<int> > & m, EHH12_ehh_data* ehhdata);

        void init_output_and_querymap(){
            this->output.resize(hm->mapData->nloci);
            //iterate over all loci
            for(int i = 0; i < hm->mapData->nloci; i++){
                this->locus_query_map[hm->mapData->mapEntries[i].locusName] = i;
                this->locus_query_map[to_string(hm->mapData->mapEntries[i].physicalPos)] = i;
            }
        }
        int queryFound(string query){
            int queryLoc = -1;
            if (locus_query_map.count(query)>0){
                queryLoc =  locus_query_map[query];
            }

            if (queryLoc < 0)
            {
                cerr << "ERROR: Could not find specific locus query, " << query << ", in data.\n";
                return -1;
            }else{
                cerr << "Found " << query << " in data.\n";
            }
            return  queryLoc;
        }


        void findMaxTwo(int* arr, int n, int &max1, int &max2) {
            max1 = 0;
            max2 = 0;
            for (int i = 0; i < n; i++) {
                if (arr[i]  > max1)
                {
                    max2 = max1;
                    max1 = arr[i] ;
                }
                else if (arr[i]  > max2)
                {
                    max2 = arr[i] ;
                }
            }
        }

        


        

};


#endif