#ifndef __SELSCAN_EHH_H__
#define __SELSCAN_EHH_H__
#include <unordered_map>
#include "selscan-stats.h"

class EHH : public SelscanStats{
    public:
        EHH(const std::unique_ptr<HapMap>&  hm, param_main& params,  ofstream* flog,  ofstream* fout) : SelscanStats(hm, params,  flog,  fout){      
        }
        void calc_single_ehh(string query);
        
    private:
        void calc_ehh_unidirection(int locus, unordered_map<unsigned int, vector<unsigned int> > & m, bool downstream);
};


#endif