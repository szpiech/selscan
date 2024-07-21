#ifndef __SELSCAN_EHH12_H__
#define __SELSCAN_EHH12_H__

#include "selscan-stats.h"

class EHH12 : public SelscanStats{
    public:
        EHH12(const std::unique_ptr<HapMap>&  hm, param_main& params,  ofstream* flog,  ofstream* fout) : SelscanStats(hm, params,  flog,  fout){      
        }
        void calc_single_ehh(string query);
        
    private:
        void calc_ehh_unidirection(int locus, unordered_map<int, vector<int> > & m, bool downstream);
};


#endif