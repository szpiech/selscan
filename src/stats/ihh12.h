#ifndef __IHH12_H__
#define __IHH12_H__

#include "selscan-stats.h"
#include <unordered_map>

using namespace std;

class IHH12_ehh_data{
    public:
    double prev_ehh12_before_norm = twice_num_pair(nhaps);
    double curr_ehh12_before_norm = 0;

    // double prev_ehh12d1_before_norm = 0;
    // double curr_ehh12d1_before_norm = 0;

    double prev_ehh_before_norm = twice_num_pair(nhaps);
    double curr_ehh_before_norm = 0;

    uint64_t n_c[2] = {0,0};
    int nhaps;
    int totgc = 0;
    int* group_count;
    bool* isDerived;
    int* group_id;

    vector<unsigned int> v;

    ~IHH12_ehh_data(){
        delete[] group_count;
        delete[] isDerived;
        delete[] group_id;
    }
    inline unsigned int square_alt(int n){
        return n*n;
    }

    inline uint64_t twice_num_pair(int n){
        // if(n < 2){
        //     return 0;
        // }
        // return 2*nCk(n, 2);
        return n*n - n;
    }
    void init(int nhaps, vector <unsigned int>& positions ){
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
        n_c[0] = nhaps - v.size();
        n_c[1] = v.size();
    }

    void initialize_core(bool ALT){
        if(n_c[1] == 0 || n_c[0] == 0){
            cout<<"WARNING: cannot use this as core"<<endl;
        }

        group_count[1] = n_c[1];
        group_count[0] = n_c[0];
        totgc+=2;
        for (int set_bit_pos : v){
            isDerived[set_bit_pos] = true;
            group_id[set_bit_pos] = 1;
        }

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
        IHH12(const std::unique_ptr<HapMap>&  hm, param_main& params,  ofstream* flog,  ofstream* fout) : SelscanStats(hm, params,  flog,  fout){}
        void main();
        void findMaxTwo(int* arr, int n, int &max1, int &max2);
        //void findMaxK(int* arr, int n, int &max1, int &max2, int k);

    private:
        static pthread_mutex_t mutex_log;
        double* ihh12;

        int max_extend;

        void calc_stat_at_locus(int locus, unordered_map<unsigned int, vector<unsigned int> >& m);
        void calc_ehh_unidirection(int locus, unordered_map<unsigned int, vector<unsigned int> > & m, bool downstream);
        void static thread_main(int tid, unordered_map<unsigned int, vector<unsigned int> >& m, IHH12* obj);

        void updateEHH_from_split(const unordered_map<unsigned int, vector<unsigned int> > & m, IHH12_ehh_data* p);
};





#endif