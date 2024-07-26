#ifndef __SELSCAN_EHH_H__
#define __SELSCAN_EHH_H__
#include <unordered_map>
#include "selscan-stats.h"

class EHHOutput{
    public:
        bool print = false;
        double ehh = 0;
        double ehh0 = 0;
        double ehh1 = 0;
        double gdist;
        int pdist;
        //for unphased
        double cehh1 = 0;
        double cehh0 = 0;
        

        // void print(){
        //     (*fout) << std::fixed <<   pdist  << "\t"
        //         <<  gdist<< "\t"
        //         << ehh1 << "\t"
        //         << ehh0 << "\t";

        //     if(downstream){
        //         (*fout) << std::fixed <<   -int(hm->mapData->mapEntries[locus].physicalPos -  hm->mapData->mapEntries[i].physicalPos)  << "\t"
        //         <<  -(hm->mapData->mapEntries[locus].geneticPos -  hm->mapData->mapEntries[i].geneticPos)<< "\t"
        //         << current_derived_ehh << "\t"
        //         << current_ancestral_ehh << "\t";
        //     }else{
        //         (*fout) << std::fixed <<   hm->mapData->mapEntries[i].physicalPos -  hm->mapData->mapEntries[locus].physicalPos  << "\t"
        //         <<  hm->mapData->mapEntries[i].geneticPos -  hm->mapData->mapEntries[locus].geneticPos<< "\t"
        //         << current_derived_ehh << "\t"
        //         << current_ancestral_ehh << "\t";
        //     }
        // }
};


class EHH : public SelscanStats{
    public:
        EHH(const std::unique_ptr<HapMap>&  hm, param_main& params,  ofstream* flog,  ofstream* fout) : SelscanStats(hm, params,  flog,  fout){    
            output.resize(hm->mapData->nloci);
        }
        void calc_single_ehh(string query);
        
    private:

        vector<EHHOutput> output;
        void calc_ehh_unidirection(int locus, unordered_map<int, vector<int> > & m, bool downstream);
        void calc_ehh_unidirection(int locus, bool downstream);

};


#endif