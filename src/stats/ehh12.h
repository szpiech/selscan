#ifndef __SELSCAN_EHH12_H__
#define __SELSCAN_EHH12_H__

#include "selscan-stats.h"

class EHH12 : public SelscanStats{
    public:
        ofstream flog_obj;
        ofstream fout_obj;
        EHH12(const std::unique_ptr<HapMap>&  hm, param_main& params) : SelscanStats(hm, params){    
            flog = &flog_obj;
            fout = &fout_obj;
            init_flog_fout("ehh12");
        }
        void calc_single_ehh(string query);
        
    private:
        void calc_ehh_unidirection(int locus, unordered_map<int, vector<int> > & m, bool downstream);
};


#endif