#ifndef __SELSCAN_IHS_H__
#define __SELSCAN_IHS_H__

#include "../selscan-maintools.h"
#include <unordered_map>

using namespace std;



class IHS : public MainTools{
    public:
   
        IHS(HapMap& hm, param_main& params,  ofstream* flog,  ofstream* fout) : MainTools(hm, params,  flog,  fout){  
        }
        void ihs_main(); //thread_ihs


    private:
        static pthread_mutex_t mutex_log;
        double* iHH0;
        double* iHH1;
        double* iHH2;

        double* ciHH0;
        double* ciHH1;
        double* ciHH2;

        //void static thread_ihs(int tid, IHS* ehh_obj);
        void static thread_ihs(int tid, unordered_map<unsigned int, vector<unsigned int> >& m, IHS* ehh_obj);

        
        void calc_ehh_unidirection_ihs(int locus, bool downstream);
        void calc_ehh_unidirection_ihs_bitset(int locus, bool downstream, unordered_map<unsigned int, vector<unsigned int> >& m);

        void calc_ehh_unidirection_ihs_unphased(int locus, bool downstream);
        double calc_ihs_unphased(int locus);

        //unphased_ihs        
        void updateEHH_from_split( map<int, vector<int> >& m, int* group_count, int* group_id, int& totgc, uint64_t* ehh_before_norm, uint64_t* cehh_before_norm, bool* is1, bool* is2);

        void calc_ihh(int locus, unordered_map<unsigned int, vector<unsigned int> >& m);  


        //helpers
        string getOrder(uint64_t n_c2, uint64_t n_c1, uint64_t n_c0);
};

#endif
