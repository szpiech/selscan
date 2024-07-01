#ifndef __XPIHH_H__
#define __XPIHH_H__

#include "../selscan-maintools.h"
#include <unordered_map>

using namespace std;

class XPIHH_ehh_data{
    public:
    uint64_t prev_ehh_before_norm = 0;
    uint64_t curr_ehh_before_norm = 0;
    uint64_t n_c[2] = {0,0};
    int nhaps;
    int totgc = 0;
    int* group_count;
    bool* isDerived;
    int* group_id;

    vector<unsigned int> v;
    vector<unsigned int> v2;
    int nhaps_p1;

    ~XPIHH_ehh_data(){
        delete[] group_count;
        delete[] isDerived;
        delete[] group_id;
    }
    inline unsigned int square_alt(int n){
        return n*n;
    }

    inline double twice_num_pair(int n){
        if(n < 2){
            return 0;
        }
        return 2*nCk(n, 2);
        //return n*n - n;
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

    void init_fix_for_pooled(vector <unsigned int>& positions, int nhaps_p1){
        v2 = positions;
        n_c[1] = v.size() + v2.size();
        n_c[0] = nhaps - n_c[1] ;   
        this->nhaps_p1 = nhaps_p1;
    }

    void initialize_core(bool ALT){
        if(n_c[1]==0){    //none set
            group_count[0] = nhaps;
            totgc+=1;
            curr_ehh_before_norm = (ALT? square_alt(n_c[0]) : twice_num_pair(n_c[0]));
        }else if (n_c[1]==nhaps){ // all set
            group_count[0] = nhaps;
            totgc+=1;
            for (int set_bit_pos : v){
                isDerived[set_bit_pos] = true;
            }
            curr_ehh_before_norm = (ALT? square_alt(n_c[1]) : twice_num_pair(n_c[1]));
        }else{
            group_count[1] = n_c[1];
            group_count[0] = n_c[0];
            totgc+=2;
            for (int set_bit_pos : v){
                isDerived[set_bit_pos] = true;
                group_id[set_bit_pos] = 1;
            }
            curr_ehh_before_norm = (ALT?  square_alt(n_c[0]) +  square_alt(n_c[1]) : twice_num_pair(n_c[0]) + twice_num_pair(n_c[1]));
        }
    }

    void initialize_core_pooled(bool ALT){
        if(n_c[1] == 0){    //none set
            group_count[0] = nhaps;
            totgc+=1;
            curr_ehh_before_norm = (ALT? square_alt(n_c[0]) : twice_num_pair(n_c[0]));
        }else if (n_c[1] == nhaps){ // all set
            group_count[0] = nhaps;
            totgc+=1;
            for (int set_bit_pos : v){
                isDerived[set_bit_pos] = true;
            }
            for (int set_bit_pos : v2){
                isDerived[set_bit_pos + nhaps_p1] = true;
            }
             curr_ehh_before_norm = (ALT? square_alt(n_c[1]) : twice_num_pair(n_c[1]));
        }else{
            group_count[1] = n_c[1];
            group_count[0] = n_c[0];
            totgc+=2;
            for (int set_bit_pos : v){
                isDerived[set_bit_pos] = true;
                group_id[set_bit_pos] = 1;
            }
            for (int set_bit_pos : v2){
                isDerived[set_bit_pos + nhaps_p1] = true;
                group_id[set_bit_pos + nhaps_p1] = 1;
            }
            curr_ehh_before_norm = (ALT?  square_alt(n_c[0]) +  square_alt(n_c[1]) : twice_num_pair(n_c[0]) + twice_num_pair(n_c[1]));
        }
    }
}; 

class XPIHH : public MainTools{
    public:
        XPIHH(HapMap& hm, param_main& params,  ofstream* flog,  ofstream* fout) : MainTools(hm, params,  flog,  fout){}
        void xpihh_main();

        

    private:
        static pthread_mutex_t mutex_log;
        double* ihh_p1;
        double* ihh_p2;

        void calc_xpihh(int locus, unordered_map<unsigned int, vector<unsigned int> >& m);
        void calc_ehh_unidirection(int locus, unordered_map<unsigned int, vector<unsigned int> > & m, bool downstream);
        void static thread_xpihh(int tid, unordered_map<unsigned int, vector<unsigned int> >& m, XPIHH* obj);

        void updateEHH_from_split(const unordered_map<unsigned int, vector<unsigned int> > & m, XPIHH_ehh_data* p);

};





#endif