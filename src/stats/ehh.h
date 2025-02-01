#ifndef __SELSCAN_EHH_H__
#define __SELSCAN_EHH_H__
#include <unordered_map>
#include "selscan-stats.h"
#include "ihs.h"

class EHHOutput{
    public:
        bool print = false;
        double ehh = 0;
        double ehh0 = 0;
        double ehh1 = 0;
        double gdist;
        int pdist;
        //for unphased
        double ehh2 = 0;
        double cehh2 = 0;
        double cehh0 = 0;
        double cehh = 0;
        
};

class EHH : public IHS{
    public:
        EHH(const std::unique_ptr<HapMap>&  hm, param_main& params) : IHS(hm, params){ }
        void main(string query);

    private:
        void calc_ehh_unidirection(int locus, bool downstream);
        void calc_ehh_unidirection_unphased(int locus, bool downstream);

        vector<EHHOutput> output;
        map<string, int> locus_query_map;

        void init_output_and_querymap(){
            this->output.resize(hm->mapData->nloci);
            for(int i = 0; i < hm->mapData->nloci; i++){
                this->locus_query_map[hm->mapData->mapEntries[i].locusName] = i;
                this->locus_query_map[to_string(hm->mapData->mapEntries[i].physicalPos)] = i;
            }
        }
        
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
                return -1;
            }else{
                cerr << "Found " << query << " in data.\n";
            }
            return  queryLoc;
        }

};


#endif