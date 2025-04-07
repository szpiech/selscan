#ifndef __SELSCAN_IHS_H__
#define __SELSCAN_IHS_H__

#include "selscan-stats.h"
#include "../thread_pool.h"
#include <unordered_map>

using namespace std;


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



    struct IhhComponents {
        double derived_right;
        double derived_left;
        double ancestral_right;
        double ancestral_left;
    };
class IHS: public SelscanStats{


    
    public:
        IHS(const std::unique_ptr<HapMap>&  hm, param_main& params) : SelscanStats(hm, params){  

        }
        void main(); 
        pair<double, double> calc_ihh1(int locus);  
        IhhComponents calc_ihh1_details(int locus);
        pair<double, double> infer_missing(int locus);  
        std::mutex mutex_log;
        //std::unique_ptr<std::mutex> mutex_log = std::make_unique<std::mutex>();


    protected:
        void updateEHH_from_split_unphased( unordered_map<int, vector<int> >& m, int* group_count, int* group_id, int& totgc, double* ehh_before_norm, double* cehh_before_norm, bool* is1, bool* is2, int* group_core);
        string getOrder(uint64_t n_c2, uint64_t n_c1, uint64_t n_c0);

    private:
        //static pthread_mutex_t mutex_log;

        

        int max_extend;

        //phased_ihs
        pair<double, double> calc_ehh_unidirection(int locus, bool downstream);

        //unphased_ihs  
        pair<double, double> calc_ehh_unidirection_unphased(int locus, bool downstream, double& cihh2, double& cihh0);

        // missing support
        pair<double, double> calc_ehh_unidirection_missing(int locus, bool downstream);
        pair<double, double> calc_ehh_unidirection_unphased_missing(int locus, bool downstream, double& cihh2, double& cihh0);
};

#endif
