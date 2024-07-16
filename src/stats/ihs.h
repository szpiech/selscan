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
            if(p.CALC_NSL){
                this->max_extend = ( p.MAX_EXTEND_NSL <= 0 ) ? physicalDistance_from_core(0,hm->hapData->nloci-1,false) : p.MAX_EXTEND_NSL;
            }else{
                this->max_extend = ( p.MAX_EXTEND <= 0 ) ? physicalDistance_from_core(0,hm->hapData->nloci-1,false) : p.MAX_EXTEND;
            }
            
        }
        void main(); //thread_ihs
        void main_old(); //thread_ihs
        //ThreadPool* pool;
        ~IHS(){
            //delete pool;
        }

        // double* iHH0;
        // double* iHH1;
        // double* iHH2;
        
        pair<double, double> calc_ehh_unidirection(int locus, bool downstream);
        pair<double, double> calc_ihh1(int locus);  

        
    private:
        static pthread_mutex_t mutex_log;
        int max_extend;

        void static thread_ihs(int tid, IHS* ehh_obj,  double& iHH1, double& iHH0);
        pair<double, double> calc_ehh_unidirection_bitset(int locus, bool downstream);

        //unphased_ihs  
        pair<double, double> calc_ehh_unidirection_unphased(int locus, bool downstream, double& cihh2, double& cihh0);
        double get_ihs_unphased(int locus);
        void updateEHH_from_split_unphased( unordered_map<int, vector<int> >& m, int* group_count, int* group_id, int& totgc, double* ehh_before_norm, double* cehh_before_norm, bool* is1, bool* is2, int* group_core);
        
        //unphased_ihs helpers
        string getOrder(uint64_t n_c2, uint64_t n_c1, uint64_t n_c0);

};



#endif
