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
        // double cehh1 = 0;
        // double cehh0 = 0;
        
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
        EHH(const std::unique_ptr<HapMap>&  hm, param_main& params) : SelscanStats(hm, params){    
            output.resize(hm->mapData->nloci);
            //iterate over all loci
            for(int i = 0; i < hm->mapData->nloci; i++){
                locus_query_map[hm->mapData->mapEntries[i].locusName] = i;
                locus_query_map[to_string(hm->mapData->mapEntries[i].physicalPos)] = i;
            }
            init_global_fout("ehh");
        }
        ~EHH(){
        }
        void calc_single_ehh(string query);
        
    private:
        map<string, int> locus_query_map;
        vector<EHHOutput> output;
        void calc_ehh_unidirection(int locus, bool downstream);

        /***
         * @param query: Locus name
         * @returns locus ( in range [0 to nloci) ) after integrity check
        */
        int queryFound(string query){
            int queryLoc = -1;
            if (locus_query_map.count(query)>0){
                queryLoc =  locus_query_map[query];
            }

            if (queryLoc < 0)
            {
                cerr << "ERROR: Could not find specific locus query, " << query << ", in data.\n";
                //cerr << "WARNING: We filtered all sites below MAF cutoff. If you wish to calculate EHH for all sites, either change --maf or set --keep-low-freq." << endl;
                return -1;
                //throw 1;
            }else{
                cerr << "Found " << query << " in data.\n";
                // double queryFreq = hapData.calcFreq(queryLoc);
                // if (hapData.SKIP && (queryFreq < hapData.MAF || 1 - queryFreq < hapData.MAF))
                // {
                //     cerr << "ERROR: EHH not calculated for " << query << ". MAF < " << hapData.MAF << ".\n";
                //     cerr << "\tIf you wish to calculate EHH at this locus either change --maf or set --keep-low-freq.\n";
                //     throw 1;
                // }
                // else if (!hapData.SKIP && (queryFreq == 0 || queryFreq == 1)){
                //     cerr << "ERROR: EHH not calculated for " << query << ". Frequency = " << queryFreq << ".\n";
                //     throw 1
                // }
                // else
                // {
                //     cerr << "Found " << query << " in data.\n";
                // }
            }
            return  queryLoc;
        }

};


#endif