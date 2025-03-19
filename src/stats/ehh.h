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
        double gdist = 0;
        int pdist = 0;
        //for unphased
        double ehh2 = 0;
        double cehh2 = 0;
        double cehh0 = 0;
        double cehh = 0;
        
};

class EHH : public SelscanStats {
    public:
        EHH(const std::unique_ptr<HapMap>&  hm, param_main& params) : SelscanStats(hm, params){ }
        void main(string query);

    private:
        void calc_ehh_unidirection(int locus, bool downstream);
        void calc_ehh_unidirection_unphased(int locus, bool downstream);

        vector<EHHOutput> output;
        map<string, int> locus_query_map;
        map<int, int> locus_query_map_position;

        void init_output_and_querymap(){
            this->output.resize(hm->mapData->nloci);
            for(int i = 0; i < hm->mapData->nloci; i++){
                this->locus_query_map[hm->mapData->mapEntries[i].locusName] = i;
                this->locus_query_map_position[hm->mapData->mapEntries[i].physicalPos] = i;
                //cout<<hm->mapData->mapEntries[i].locusName<<" "<<hm->mapData->mapEntries[i].physicalPos<<endl;
            }
        }

        string getOrder(uint64_t n_c2, uint64_t n_c1, uint64_t n_c0){
            string order_str;
            if(n_c2 >= n_c1 and n_c1 >= n_c0){
                order_str = "210";
            }else if(n_c1 >= n_c2 and n_c2 >= n_c0){
                order_str = "120";
            }else if(n_c2 >= n_c0 and n_c0 >= n_c1){
                order_str = "201";
            }else if(n_c1 >= n_c0 and n_c0 >= n_c2){
                order_str = "102";
            }else if(n_c0 >= n_c2 and n_c2 >= n_c1){
                order_str = "021";
            }else if(n_c0 >= n_c1 and n_c1 >= n_c2){
                order_str = "012";
            }
            return order_str;
        }
        
        void updateEHH_from_split_unphased( unordered_map<int, vector<int> >& m, int* group_count, int* group_id, int& totgc, double* ehh_before_norm, double* cehh_before_norm, bool* is1, bool* is2, int* group_core){
            for (const auto &ele : m) {
                int old_group_id = ele.first;
                int newgroup_size = ele.second.size() ;
        
                if(group_count[old_group_id] == newgroup_size || newgroup_size == 0){
                    continue;
                }
                for(int v: ele.second){
                    group_id[v] = totgc;
                }
                 
                group_count[old_group_id] -= newgroup_size;
                group_count[totgc] += newgroup_size;
                
                totgc+=1;
                
                //TODO
                if(is1[ele.second[0]]) // if the core locus for this chr has 1 (derived), then update ehh1, otherwise ehh0
                {
                    group_core[totgc-1] = 1; // cannot directly use del_update to update for unphased
                }else if(is2[ele.second[0]]) {
                    group_core[totgc-1] = 2;
                }else{
                    group_core[totgc-1] = 0;
                }   
            }
        }

        /***
         * @param query: Locus name
         * @returns locus ( in range [0 to nloci) ) after integrity check
        */
        int queryFound(string query){
            int queryLoc = -1;

            if (locus_query_map.count(query) < 1)
            {
                std::cerr << "Could not find specific locus query ID " << query << " in filtered data. Searching for POS "<< query << "... \n";   
            }else{
                queryLoc =  locus_query_map[query];
            }

            if(queryLoc!=-1){
                cerr << "Found ID " << query << " in data. Computing EHH...\n";
                return queryLoc;
            }
            
            try {
                int query_pos = std::stoi(query);  // Attempt conversion
                if (locus_query_map_position.count(query_pos) < 1){
                    throw std::exception();
                }else{
                    queryLoc =  locus_query_map_position[query_pos];
                }
            } catch (const std::exception & e) {
                std::cerr << "ERROR: Could not find specific locus query ID or POS "<< query<< " in filtered data." << std::endl;
                std::cerr <<  "If the site is filtered out because of MAF cutoff, either reduce --maf or set --keep-low-freq." << std::endl;
                exit(EXIT_FAILURE);
            } 
            
            
            if(queryLoc!=-1){
                cerr << "Found POS " << query << " in data. Computing EHH...\n";
                return queryLoc;
            }
            
                // cerr << "Available locus queries are: ";
                // //iterate over locus_query_map and print
                // for(auto const& x : locus_query_map){
                //     cerr << x.first  // string (key)
                //               << ',' ;
                // }
                // cerr<<endl;
        }



};


#endif