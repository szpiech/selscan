#ifndef __SELSCAN_EHH12_H__
#define __SELSCAN_EHH12_H__

#include "../selscan-maintools.h"

class EHH12 : public MainTools{
    public:
        EHH12(HapMap& hm, param_main& params,  ofstream* flog,  ofstream* fout) : MainTools(hm, params,  flog,  fout){      
        }
        void calc_single_ehh(string query);
        
    private:
        void calc_ehh_unidirection(int locus, unordered_map<unsigned int, vector<unsigned int> > & m, bool downstream);
};


#endif