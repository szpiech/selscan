#ifndef __SELSCAN_IHS_H__
#define __SELSCAN_IHS_H__

#include "selscan-stats.h"
#include "../thread_pool.h"
#include <unordered_map>

using namespace std;

class IHS: public SelscanStats{
    public:
        IHS(const std::unique_ptr<HapMap>&  hm, param_main& params,  ofstream* flog,  ofstream* fout) : SelscanStats(hm, params,  flog,  fout){  
            //pool = new ThreadPool(numThreads);
        }
        void main(); //thread_ihs
        void main_old(); //thread_ihs
        //ThreadPool* pool;
        ~IHS(){
            //delete pool;
        }

        double* iHH0;
        double* iHH1;
        double* iHH2;
        
        pair<double, double> calc_ehh_unidirection(int locus, bool downstream);
        pair<double, double> calc_ihh1(int locus);  
        pair<double, double>  calc_ihh(int locus);  

        
    private:
        static pthread_mutex_t mutex_log;
      
        double* ciHH0;
        double* ciHH1;
        double* ciHH2;

        void static thread_ihs(int tid, IHS* ehh_obj);
        //pair<double, double>  static  thread_ihs(int tid, unordered_map<unsigned int, vector<unsigned int> >& m, IHS* ehh_obj);
        
        pair<double, double> calc_ehh_unidirection_bitset(int locus, bool downstream);

        //unphased_ihs  
        void calc_ehh_unidirection_unphased(int locus, bool downstream);
        double get_ihs_unphased(int locus);
        void updateEHH_from_split_unphased( map<int, vector<int> >& m, int* group_count, int* group_id, int& totgc, uint64_t* ehh_before_norm, uint64_t* cehh_before_norm, bool* is1, bool* is2);
        //unphased_ihs helpers
        string getOrder(uint64_t n_c2, uint64_t n_c1, uint64_t n_c0);


        void getNextSetBits(vector<unsigned int>*, bool, int);

};



#endif
