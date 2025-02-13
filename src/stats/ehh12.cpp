#include "ehh12.h"




/**
 * Calculate EHH in only one direction until cutoff is hit - upstream or downstream
*/
void EHH12::calc_ehh12_unidirection(int locus, bool downstream){

    double ihh12 = 0;
    unordered_map<int, vector<int> > m;

    EHH12_ehh_data ehhdata;
    int numSnps = hm->hapData->nloci; // must be same for both hapData and hapData2

    ehhdata.init(hm->hapData->nhaps, get_all_1s(locus));
    ehhdata.initialize_core(p.ALT);


    output[locus].ehh1 = ehhdata.curr_ehh_before_norm / twice_num_pair(ehhdata.nhaps);
    output[locus].ehh12 = ehhdata.curr_ehh12_before_norm  / twice_num_pair(ehhdata.nhaps);
    output[locus].ehh2d1 = ehhdata.curr_ehh12d1_before_norm / twice_num_pair(ehhdata.nhaps);

    int i = locus;
    while(true){

        if(physicalDistance_from_core(i, locus, downstream) >= p.QWIN){ //check if currentLocus is beyond 1Mb
            break;
        }

        double breaking_ehh = (p.ALT ? ehhdata.curr_ehh_before_norm / square_alt(ehhdata.nhaps) : ehhdata.curr_ehh_before_norm /twice_num_pair(ehhdata.nhaps));
        if(breaking_ehh <= p.EHH_CUTOFF){
            /*DEBUG :
            if(downstream){
                cout<<"breaking ehh down: "<<locus<<" "<<breaking_ehh<<" "<<ihh12[locus]<<endl;
            }else{
                cout<<"breaking ehh up: "<<locus<<" "<<breaking_ehh<<" "<<ihh12[locus]<<endl;
            }
             */
            break;
        }

        bool breakReachedEdge = false;
        breakReachedEdge = downstream? (i == 0) : (i == numSnps-1);
        if(breakReachedEdge){ //nextLocus < 0 || nextLocus >= numSnps to avoid going to negative
            (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << p.EHH_CUTOFF << ". " << endl;
            break;
        }

        //at this point we ensured that nextLocus (i-1 or i+1) is within bounds
        i = downstream? i-1 : i+1; //proceed to next locus
        // int physicalDistance_old = physicalDistance(i,downstream); //double check if i or nextlocus
        // if (physicalDistance_old > p.MAX_GAP)
        // {
        //     (*flog) << "WARNING: Reached a gap of " << physicalDistance_old << "bp > " << p.MAX_GAP << "bp. Skipping calculation at position " <<  hm->mapData->mapEntries[i].physicalPos << " id: " <<  hm->mapData->mapEntries[i].locusName << "\n";
        //     return;
        // }

        ehhdata.v = get_all_1s(i);
        ACTION_ON_ALL_SET_BITS(ehhdata.v, {
            int old_group_id = ehhdata.group_id[set_bit_pos];
            m[old_group_id].push_back(set_bit_pos); 
        });
        
        updateEHH_from_split(m, &ehhdata);
        m.clear();

        //if(!downstream) cout<<current_ehh1<<" "<<previous_ehh1<<" "<<current_ehh12<<" "<<previous_ehh12<<" "<<ihh12[locus] <<" "<<scale * distance<< endl;

        output[i].ehh1 = ehhdata.curr_ehh_before_norm / twice_num_pair(ehhdata.nhaps);
        output[i].ehh12 = ehhdata.curr_ehh12_before_norm  / twice_num_pair(ehhdata.nhaps);
        output[i].ehh2d1 = ehhdata.curr_ehh12d1_before_norm / twice_num_pair(ehhdata.nhaps);

        if (physicalDistance_from_core(i, locus, downstream) >= p.QWIN) break;
        
    }
}

void EHH12::main(string query){
    init_global_fout("ehh12");
    init_output_and_querymap();
    
    int numSnps = hm->mapData->nloci;
    int numHaps = hm->hapData->nhaps;

    int locus = this->queryFound(query);
    if(locus == -1){
        return;
    }

    if(p.UNPHASED){
        throw "ERROR: --unphased and --ehh12 not compatible";
        exit(EXIT_FAILURE);
    }else{
        calc_ehh12_unidirection(locus, true); // downstream
        calc_ehh12_unidirection(locus, false); // upstream

        (*fout) << std::fixed <<   "pdist\tgdist\tehh1\tehh12\tehh2d1" << "\n";
        for (int i = 0; i < numSnps; i++){
            if(!output[i].print){
                continue;
            }

            fout->precision(6);
            (*fout) << std::fixed <<   output[i].pdist  << "\t"
            <<  output[i].gdist << "\t"
            << output[i].ehh1 << "\t"
            << output[i].ehh12  << "\t"
            << output[i].ehh2d1  << "";
            (*fout) << endl;
        }
    }
}



void EHH12::updateEHH_from_split(const unordered_map<int, vector<int> > & m, EHH12_ehh_data* ehhdata){
    double sum_del_update = 0;
    for (const auto &ele : m) {
        int old_group_id = ele.first;
        int newgroup_size = ele.second.size() ;

        if(ehhdata->group_count[old_group_id] == newgroup_size || newgroup_size == 0){ //if a group becomes empty, we don't increment num_groups, we just reuse that group
            continue;
        }

        for(int v: ele.second){
            ehhdata->group_id[v] = ehhdata->totgc;
        }

        double del_update = -twice_num_pair(ehhdata->group_count[old_group_id]) + twice_num_pair(newgroup_size) + twice_num_pair(ehhdata->group_count[old_group_id] - newgroup_size);
        if(p.ALT){
            del_update = -square_alt(ehhdata->group_count[old_group_id]) +   square_alt(newgroup_size) + square_alt(ehhdata->group_count[old_group_id] - newgroup_size);
        }

        ehhdata->group_count[old_group_id] -= newgroup_size;
        ehhdata->group_count[ehhdata->totgc] += newgroup_size;
        ehhdata->totgc += 1;

        sum_del_update += del_update;
    }
    ehhdata->curr_ehh_before_norm += sum_del_update;
    int top1, top2;
    findMaxTwo(ehhdata->group_count, ehhdata->totgc, top1, top2);

    double firstFreq = (top1 > 1) ? twice_num_pair(top1) : 0;
    double secondFreq =(top2 > 1) ?  twice_num_pair(top2): 0;
    double comboFreq = ((top1 + top2) > 1) ? twice_num_pair((top1 + top2)) : 0;
    double normfac = twice_num_pair(ehhdata->nhaps);

    ehhdata->curr_ehh12_before_norm = ehhdata->curr_ehh_before_norm  - firstFreq - secondFreq + comboFreq;

    ehhdata->curr_ehh12d1_before_norm = ehhdata->curr_ehh_before_norm  - firstFreq;

    //ehhdata->curr_ehh12_before_norm *= normfac;
    //cout<<"t1 t2 "<<top1<<" "<<top2<<" "<<firstFreq/normfac<<" "<<secondFreq/normfac<<" "<< comboFreq/normfac<<endl;
}