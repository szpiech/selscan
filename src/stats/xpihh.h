#ifndef __XPIHH_H__
#define __XPIHH_H__

#include "selscan-stats.h"
#include <unordered_map>

using namespace std;

struct XPIHH_ehh_data{
    double prev_ehh_before_norm = 0;
    double curr_ehh_before_norm = 0;
    int n_c[2] = {0,0};

    int nhaps;
    int totgc = 0;
    int* group_count;
    bool* isDerived;
    int* group_id;
    double normalizer;
    bool LOW_MEM;


    ~XPIHH_ehh_data();
    // inline int square_alt(int n);
    // inline double twice_num_pair(int n);
    // void init(int nhaps, const vector <int>* positions, bool ALT, bool WAGH );
    // void init_for_pooled(const vector <int>* positions,const vector <int>* positions2, int nhaps_p1, bool ALT, bool WAGH );
    // void initialize_core(const vector <int>* v);
    // void initialize_core_pooled(const vector <int>* v, const vector <int>* v2,  int nhaps_p1); 
}; 


class XPIHH : public SelscanStats{
    public:
        XPIHH(const std::unique_ptr<HapMap>&  hm, param_main& params, ofstream* flog, ofstream* fout) : SelscanStats(hm, params,  flog,  fout){
            if(p.CALC_XPNSL){
                this->max_extend = ( p.MAX_EXTEND_NSL <= 0 ) ? physicalDistance_from_core(0,hm->hapData->nloci-1,true) : p.MAX_EXTEND_NSL;
            }else{
                this->max_extend = ( p.MAX_EXTEND <= 0 ) ? physicalDistance_from_core(0,hm->hapData->nloci-1,true) : p.MAX_EXTEND;
            }
        }
        void main();

    private:
        static pthread_mutex_t mutex_log;
        int max_extend;
        pair<double, double> calc_xpihh(int locus);
        pair<double, double> calc_ehh_unidirection(int locus, bool downstream);
        void updateEHH_from_split(const unordered_map<int, vector<int> > & m, XPIHH_ehh_data* p);
};


#endif