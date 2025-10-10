#include "ihs.h"
#include <iomanip> 
#include <set>
#include <algorithm>
#include <cassert>
#include <pthread.h>

/**
 * Update extended haplotype homozygosity (EHH) groups after splitting chromosomes (unphased).
 *
 * @param m                  A map from old group ID to a list of chromosome indices that are being split into a new group.
 * @param group_count        Array of counts of individuals in each group.
 * @param group_id           Array mapping each chromosome to its current group ID.
 * @param totgc              Total number of unique groups created so far. Will be incremented.
 * @param ehh_before_norm    Placeholder for EHH values before normalization-by-tot-core-pairs (not updated here).
 * @param cehh_before_norm   Placeholder for complement EHH values before normalization-by-tot-core-pairs (not updated here).
 * @param is1                Boolean array indicating if each chromosome has the derived allele at the core SNP.
 * @param is2                Boolean array indicating if each chromosome has the ancestral allele at the core SNP.
 * @param group_core         Array tracking the core allele (1 = derived, 2 = ancestral, 0 = unknown) for each group.
 */
void IHS::updateEHH_from_split_unphased( unordered_map<int, vector<int> >& m, int* group_count, int* group_id, int& totgc, bool* is1, bool* is2, int* group_core){
    for (const auto &ele : m) {
        int old_group_id = ele.first;
        int newgroup_size = ele.second.size() ;

        // Skip empty splits or unchanged groups
        if(group_count[old_group_id] == newgroup_size || newgroup_size == 0){
            continue;
        }

        // Assign new group ID to each chromosome in the split
        for(int v: ele.second){
            group_id[v] = totgc;
        }
         
        // Update group counts: reduce old group size, increase new group size
        group_count[old_group_id] -= newgroup_size;
        group_count[totgc] += newgroup_size;
        
        // Increment total group counter for the next group
        totgc+=1;
        
        // Assign the core allele state for the new group
        if(is1[ele.second[0]]) 
        {
            group_core[totgc-1] = 1; // heterozygous
        }else if(is2[ele.second[0]]) {
            group_core[totgc-1] = 2; // homozygous derived
        }else{
            group_core[totgc-1] = 0; // homozygous ancestral
        }   
    }
}


/**
 * Compute ciHH2, ciHH0, iHH2, iHH0 in one direction from a given core SNP.
 *
 * This function calculates the unphased extended haplotype homozygosity (EHH)
 * in one direction (either upstream or downstream) starting from a core locus.
 * It accumulates the integrated haplotype homozygosity (iHH) for chromosomes
 * carrying the derived (core = 1) and ancestral (core = 0) alleles separately,
 * and returns a pair of iHH values (derived, ancestral).
 *
 * @param locus      Index of the core SNP from which to begin EHH calculation.
 * @param downstream Direction of integration: true for downstream (right), false for upstream (left).
 * @param cihh2      Used for additional Output: ciHH for derived homozygous.
 * @param cihh0      Used for additional Output: ciHH for ancestral homozygous.
 * @return           A pair: {iHH_derived, iHH_ancestral}
 */
OutputUnphasedIHH IHS::calc_ehh_unidirection_unphased(int locus, bool downstream){
    OutputUnphasedIHH skipLocusOutput = {SKIP_LOCUS_VALUE, SKIP_LOCUS_VALUE, SKIP_LOCUS_VALUE, SKIP_LOCUS_VALUE};
    std::unique_ptr<std::unordered_map<int, std::vector<int>>> mp(new std::unordered_map<int, std::vector<int>>());
    unordered_map<int, vector<int> >& m = (* mp);

    int numHaps = hm->hapData->nhaps;

    //double prev_ehh_before_norm[3];
    double curr_ehh_before_norm[3];
    //double  prev_cehh_before_norm[3]; 
    double  curr_cehh_before_norm[3];

    double prev_ehh[3] = {1, 1, 1};
    double prev_cehh[3] = {1, 1, 1};


    int n_c[3] = {0,0,0};

    int* group_count = new int[numHaps];
    int* group_id =  new int[numHaps];
    int* group_core =  new int[numHaps];
    bool* is1 =  new bool[numHaps];
    bool* is2 =  new bool[numHaps];

    double iHH[3] = {0,0,0};
    double ciHH[3] = {0,0,0};
    
    for(int i = 0; i<numHaps; i++){ //hopefully this will be vectorized with compile time flags
        group_count[i] = 0;
        is1[i] = false;
        is2[i] = false;
    }

    int totgc=0;

    n_c[1] = hm->hapData->get_n_c1(locus); //hm->hapData->hapEntries[locus].hapbitset->num_1s;
    n_c[2] = hm->hapData->get_n_c2(locus); // hm->hapData->hapEntries[locus].xorbitset->num_1s;
    n_c[0] = numHaps - n_c[1] - n_c[2]; // assume no missing
    

    // if(hm->mapData->mapEntries[locus].physicalPos==2455){
    //     cout<<"n_c1: "<<n_c[1]<<" n_c2: "<<n_c[2]<<" n_c0: "<<n_c[0]<<endl;
    //     //HANDLE_ERROR("");
    // }
    if(n_c[1] + n_c[2] + n_c[0] != numHaps){
        HANDLE_ERROR("n_c1 + n_c2 + n_c0 != numHaps");
    }

    string orderStr = getOrder(n_c[2], n_c[1], n_c[0]);

    double normalizer[3];
    double normalizer_not[3];
    for(int i = 0; i<3; i++){
        normalizer[i] = twice_num_pair_or_square(n_c[i], p.ALT);
        normalizer_not[i] = twice_num_pair_or_square(numHaps - n_c[i], p.ALT); 
    }

    for(int i = 0; i<3; i++){
        curr_ehh_before_norm[i] = normalizer[i];
        curr_cehh_before_norm[i] = normalizer_not[i];
        group_count[i] = n_c[orderStr[i]-'0'];
        group_core[i] = orderStr[i]-'0';
    }

    //group_count
    //[0] = most occurring
    //[1] = second most occurring
    //[2] = least occurring

    int pos_of_012[3] = {0,0,0};
    for(int i = 0; i<3; i++){
        pos_of_012[orderStr[i]-'0'] = i ;
    }

    for (int i = 0; i<numHaps; i++){
        group_id[i] = pos_of_012[0];
    }

    ACTION_ON_ALL_SET_BITS(hm->hapData->hapEntries[locus].hapbitset, {
        is1[set_bit_pos] = true;
        group_id[set_bit_pos] = pos_of_012[1];   
    });

    ACTION_ON_ALL_SET_BITS(hm->hapData->hapEntries[locus].xorbitset, {
        is2[set_bit_pos] = true;
        group_id[set_bit_pos] = pos_of_012[2]; 
    });
   
    if(group_count[0] == numHaps){ //monomorphic site
        totgc+=1;
    }else if(group_count[2] == 0 ){ //second clause is redundant
        if(group_count[0] + group_count[1] != numHaps){
            HANDLE_ERROR("gc2==0 at locus");
        }
        totgc+=2;
    }else{
        totgc+=3;
    }

    // for (int i : {0, 2}) {
    //     //we dont need i = 1
    //     prev_ehh_before_norm[i] = curr_ehh_before_norm[i];
    //     prev_cehh_before_norm[i] =  curr_cehh_before_norm[i];
    // }

    if(n_c[1] == numHaps || n_c[0] == numHaps || n_c[2] == numHaps){ 
        {
            std::lock_guard<std::mutex> lock(mutex_log);
            (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
            << " (pos " << hm->mapData->mapEntries[locus].physicalPos << ") is monomorphic. Skipping calculation at this locus. "
            << "het: " <<  n_c[1] << " hom0: " << n_c[0] << " hom1: " << n_c[2] << ".\n";
        }
        return skipLocusOutput;
    }

    // if(n_c[2] == 1){
    //     {std::lock_guard<std::mutex> lock(mutex_log);
    //     (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
    //             << " (pos " << hm->mapData->mapEntries[locus].physicalPos << ") has het count = 1. Skipping calculation at this locus.\n";
    //     }//unlock
    //     return skipLocusOutput;
    // }

    // if(n_c[0] == 1){
    //     {std::lock_guard<std::mutex> lock(mutex_log);
    //     (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
    //             << " (pos " << hm->mapData->mapEntries[locus].physicalPos << ") has ancestral homozygotes count = 1. Skipping calculation at this locus.\n";
    //     }//unlock
    //     return skipLocusOutput;
    // }

    // if(n_c[2] == 1){
    //     {std::lock_guard<std::mutex> lock(mutex_log);
    //     (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
    //             << " (pos " << hm->mapData->mapEntries[locus].physicalPos << ") has derived homozygotes count = 1. Skipping calculation at this locus.\n";
    //     }//unlock
    //     return skipLocusOutput;
    // }


    double freqHetGT = n_c[1]*1.0/numHaps;
    if (  freqHetGT > 1-p.MAF ) 
    {
        {
            std::lock_guard<std::mutex> lock(mutex_log);
            (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
                << " (pos " << hm->mapData->mapEntries[locus].physicalPos << ") has too many hets. Skipping calculation at this locus. "
                << "het: " <<  n_c[1] << " hom0: " << n_c[0] << " hom1: " << n_c[2] << ".\n";
        }//unlock
        return skipLocusOutput;
    }

    int i = locus;  // locus == core_locus
    int prev_index = locus;
    int skipped_due_to_multimaf = 0;
    while(true){ // Upstream: for ( int i = locus+1; i<all_positions.size(); i++ )
    
        if(p.CALC_IHS && !p.CALC_NSL){
            if(curr_ehh_before_norm[2]*1.0/normalizer[2] <= p.EHH_CUTOFF and curr_ehh_before_norm[0]*1.0/normalizer[0]  <= p.EHH_CUTOFF){   // or cutoff, change for benchmarking against hapbin
                DBG("Break reason for locus "<<locus<<":: EHH_CUTOFF.");
                break;
            }
        }

        bool edgeBreak = false;
        edgeBreak = nextLocOutOfBounds(i, downstream);
        if(edgeBreak) {
            {
                std::lock_guard<std::mutex> lock(mutex_log);
                (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << p.EHH_CUTOFF
                    << ". position: "<< hm->mapData->mapEntries[locus].physicalPos << " id: " << hm->mapData->mapEntries[locus].locusName ;
            }
            
            if (!p.TRUNC){
                {
                    std::lock_guard<std::mutex> lock(mutex_log);
                     (*flog) << "Skipping calculation." << endl;
                }
                return skipLocusOutput;
            }else{
                std::lock_guard<std::mutex> lock(mutex_log);
                (*flog) << endl;
            }
            break;
        }
        
        i = (downstream) ? i-1 : i+1;
        if(p.MULTI_MAF){
            if(hm->hapData->get_maf(i) < p.MAF){
                skipped_due_to_multimaf++;
                continue; // skip this locus for integration
            }
        }
        
        double distance =  geneticDistance(i, prev_index, downstream);
        if(p.CALC_NSL && !p.CALC_IHS){
            //distance = snpDistance(i, prev_index, downstream); // distance = 1; if no site was filtered out
            distance = 1;
        }
        double scale = double(p.SCALE_PARAMETER) / double(physicalDistance(i, prev_index, downstream) );
        if (scale > 1) scale = 1;
        distance *= scale;
        if (physicalDistance(i, prev_index, downstream) > p.MAX_GAP)
        {
            {
                std::lock_guard<std::mutex> lock(mutex_log);
                (*flog) << "WARNING: Reached a gap of " << physicalDistance(i, prev_index, downstream)
                        << "bp > " << p.MAX_GAP << "bp. Skipping calculation at position " << hm->mapData->mapEntries[locus].physicalPos << " id: " << hm->mapData->mapEntries[locus].locusName << "\n";
            }
            return skipLocusOutput;
        }
        prev_index = i;  // save last index that was not skipped


        /*
        if( at i it is monomorphic){
            // {std::lock_guard<std::mutex> lock(mutex_log);
            // (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
            //         << " (number " << locus << ") is monomorphic. Skipping calculation at this locus. "
            //         << "het: " <<  n_c[1] << " hom0: " << n_c[0] << " hom1: " << n_c[2] << ".\n";
            // }//unlock
            // hm->mapData->mapEntries[locus].skipLocus = true;
            // break;
            // if you wish to continue
            for(int i : {0, 2}){
                if(normalizer[i]!=0){
                    iHH[i] += (curr_ehh_before_norm[i] * 1.0 / normalizer[i]  + prev_ehh_before_norm[i] * 1.0  / normalizer[i] ) * 0.5 * distance;
                }
                if(normalizer_not[i]!=0){
                    ciHH[i] += (curr_cehh_before_norm[i] * 1.0  / normalizer_not[i] + prev_cehh_before_norm[i] * 1.0  / normalizer_not[i]) * 0.5  *  distance;
                }
            }
            continue; 
        }
        */

       if(hm->hapData->get_n_c0(i) == numHaps or hm->hapData->get_n_c1(i) == numHaps or hm->hapData->get_n_c2(i) == numHaps){ // monomorphic check, compute separately          
            DBG("WARNING: Monomorphic site at locus "<<i);
            //no need to do grouping logic  unnecessarily
            for(int i : {0, 2}){
                if(normalizer[i]!=0){
                    iHH[i] += (curr_ehh_before_norm[i] * 1.0 / normalizer[i]  + prev_ehh[i] ) * 0.5 * distance;
                }else{
                    iHH[i] += prev_ehh[i] * 0.5 * distance;
                    prev_ehh[i] = 0; // AC = 0 or 1, this ensures only first two points are used for area calculation
                }

                if(normalizer_not[i]!=0){
                    ciHH[i] += (curr_cehh_before_norm[i] * 1.0  / normalizer_not[i] + prev_cehh[i]) * 0.5  *  distance;
                }else{
                    ciHH[i] += prev_cehh[i] * 0.5 * distance;
                    prev_ehh[i] = 0; // AC = 0 or 1, this ensures only first two points are used for area calculation
                }
            }
            continue; 
        }
        
        ACTION_ON_ALL_SET_BITS(hm->hapData->hapEntries[i].xorbitset, {
            int old_group_id = group_id[set_bit_pos];
            m[old_group_id].push_back(set_bit_pos);
        });
        updateEHH_from_split_unphased(m, group_count, group_id, totgc, is1, is2, group_core);
        m.clear();

        ACTION_ON_ALL_SET_BITS(hm->hapData->hapEntries[i].hapbitset, {
            int old_group_id = group_id[set_bit_pos];
            m[old_group_id].push_back(set_bit_pos);
        });
        updateEHH_from_split_unphased(m, group_count, group_id, totgc, is1, is2, group_core);
        m.clear();

        // equivalent to calcHomozoygosity (without the normalization-by-tot-core-pairs)
        double cehh2_before_norm = 0;
        double cehh0_before_norm = 0;
        double ehh2_before_norm = 0;
        double ehh0_before_norm = 0;
        for(int x = 0; x<totgc; x++){
            long double gcsquare = twice_num_pair_or_square(group_count[x],p.ALT);
            if(group_core[x]==0){
                ehh0_before_norm += gcsquare;
                cehh2_before_norm += gcsquare;
            }else if(group_core[x]==1){
                cehh2_before_norm += gcsquare;
                cehh0_before_norm += gcsquare;
            }else{
                ehh2_before_norm += gcsquare;
                cehh0_before_norm += gcsquare;
            }
        }

        // area update for iHH0
        if(prev_ehh[0] > p.EHH_CUTOFF || (p.CALC_NSL && !p.CALC_IHS) ){ // bug fix for low MAC
            curr_ehh_before_norm[0] = ehh0_before_norm;
            if( normalizer[0]!=0 ){ 
                iHH[0] += (curr_ehh_before_norm[0]/ normalizer[0] + prev_ehh[0]) * 0.5 * distance;
                prev_ehh[0] = curr_ehh_before_norm[0] / normalizer[0];
            }else{
                iHH[0] += prev_ehh[0] * 0.5  * distance;
                prev_ehh[0] = 0; // AC = 0 or 1, this ensures only first two points are used for area calculation
            }
        }

        // area update for iHH2
        if(prev_ehh[2] > p.EHH_CUTOFF || (p.CALC_NSL && !p.CALC_IHS)){  // bug fix for low MAC
            curr_ehh_before_norm[2] = ehh2_before_norm;
            if( normalizer[2]!=0){
                iHH[2] += (curr_ehh_before_norm[2] / normalizer[2] + prev_ehh[2]) * 0.5 * distance;
                prev_ehh[2] = curr_ehh_before_norm[2]  / normalizer[2];
            }else{
                iHH[2] += prev_ehh[2] * 0.5  * distance;
                prev_ehh[2] = 0; // if normalizer is zero, we set prev_ehh to zero
            }
        }
            
        // area update for cIHH 
        if(true){  //this is how it's in selscan, does not depend on cutoff
            curr_cehh_before_norm[0] = cehh0_before_norm;
            if(normalizer_not[0]!=0){
                ciHH[0] += (curr_cehh_before_norm[0] / normalizer_not[0]  + prev_cehh[0]) * 0.5 * distance ;
                prev_cehh[0] = curr_cehh_before_norm[0] /  normalizer_not[0];
            }else{
                ciHH[0] +=  prev_cehh[0] * 0.5 * distance ;
                prev_cehh[0] = 0; // AC = 0 or 1, this ensures only first two points are used for area calculation
            }

            curr_cehh_before_norm[2] = cehh2_before_norm;
            if(normalizer_not[2]!=0){
                ciHH[2] += (curr_cehh_before_norm[2] / normalizer_not[2] + prev_cehh[2]) * 0.5 * distance ;
                prev_cehh[2] = curr_cehh_before_norm[2] / normalizer_not[2];
            }else{
                ciHH[2] +=  prev_cehh[2] * 0.5 * distance ;
                prev_cehh[2] = 0;
            }
        }

        // //If locus is monomorphic, shoot a warning and skip locus
        // //This probably isn't necessary any more
        // if ( !unphased && (numDerived == 0 || numAncestral == 0) ) 
        // {
        //     {std::lock_guard<std::mutex> lock(mutex_log);
        //     (*flog) << "WARNING: locus " << locusName[locus]
        //             << " (number " << locus + 1 << ") is monomorphic. Skipping calculation at this locus.\n";
        //     }//unlock
        //     skipLocus = 1;
        //     break;
        // }

        if(totgc == numHaps) {
            DBG("Break reason for locus "<<locus<<":: ALL_UNIQUE.");
            break;
        }

        //@DEBUG_BLOCK
        // if(locus==1){ // debugging specific locus
        //     if(downstream){
        //         cout<<locus<<":::l "<<i << " "<<curr_cehh_before_norm[2]/normalizer_not[2]<<" "<<curr_ehh_before_norm[2]/normalizer[2]<<" "<<ciHH[2]<<" "<<iHH[2]<<endl;
        //     }else{
        //         cout<<locus<<":::r "<<i << " "<<curr_cehh_before_norm[2]/normalizer_not[2]<<" "<<curr_ehh_before_norm[2]/normalizer[2]<<" "<<ciHH[2]<<" "<<iHH[2]<<endl;
        //     }
        // }

        if ((!p.CALC_NSL && p.CALC_IHS) && physicalDistance(i, locus, downstream) >= max_extend) { // max-exted is --max-extend value
            DBG("Break reason for locus "<<locus<<":: MAX_EXTEND.");
            break;
        }
        if ((p.CALC_NSL && !p.CALC_IHS) && abs(i-locus)-skipped_due_to_multimaf >= max_extend) { // max-exted is --max-extend-nsl value
            DBG("Break reason for locus "<<locus<<":: MAX_EXTEND_NSL." << i << " "<<locus<< " " << skipped_due_to_multimaf<< " " <<abs(i-locus)-skipped_due_to_multimaf << " >= " << max_extend);
            break;
        }

        // @DEBUG_BLOCK
        // if(downstream){
        //     cout<<locus<<":::l "<<i << " "<<curr_cehh_before_norm[0]/normalizer_not[0]<<" "<<curr_ehh_before_norm[0]/normalizer[0]<<" "<<ciHH[0]<<" "<<iHH[0]<<endl;
        // }else{
        //     cout<<locus<<":::r "<<i << " "<<curr_cehh_before_norm[0]/normalizer_not[0]<<" "<<curr_ehh_before_norm[0]/normalizer[0]<<" "<<ciHH[0]<<" "<<iHH[0]<<endl;
        // }
    }
    delete[] group_count;
    delete[] group_id;
    delete[] is1;
    delete[] is2;
    delete[] group_core;

    // cihh2 = ciHH[2];
    // cihh0 = ciHH[0];
    // return make_pair(iHH[2], iHH[0]);
    return {iHH[2], iHH[0],  ciHH[2],  ciHH[0]};
}


/**
 * Determine the relative order of group sizes (core2, core1, core0).
 *
 * Given counts of chromosomes in three allele groups (typically corresponding to
 * genotype states or core haplotype groups: core2, core1, and core0), this function
 * returns a string representing the order of their sizes in descending order.
 *
 * For example, if core2 > core1 > core0, the function returns "210", meaning:
 * - '2' (core2) is the largest,
 * - '1' (core1) is second,
 * - '0' (core0) is the smallest.
 *
 * @param n_c2  Count of chromosomes in group 2 (typically core2, e.g., derived).
 * @param n_c1  Count of chromosomes in group 1 (e.g., heterozygous or mixed group).
 * @param n_c0  Count of chromosomes in group 0 (typically core0, e.g., ancestral).
 * @return      A string of digits (permutation of "012") indicating the sorted order
 *              of the group sizes from largest to smallest.
 */
string IHS::getOrder(uint64_t n_c2, uint64_t n_c1, uint64_t n_c0){
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

/**
 * Calculate EHH in only one direction until cutoff is hit - upstream or downstream
*/
pair<double, double> IHS::calc_ehh_unidirection(int locus, bool downstream){

    bool OUT_EHHS = false;
    std::ofstream ehhs_out;
    string stat = p.CALC_NSL ? "nsl" : "ihs";
    string ehh_filename = p.outFilename + + "." + stat + ".ehh" + "." + std::to_string(hm->mapData->mapEntries[locus].physicalPos) + ".out";

    if(p.CALC_EHHS){
        string input = p.EHHS_RANGES; //string input = "1000-2000,5000-6000,1250000-1300000,2000000-2500000";  // to test, order doesn't matter
        vector<pair<int, int> > ranges = parse_ranges(input);
        OUT_EHHS = (is_position_in_ranges(hm->mapData->mapEntries[locus].physicalPos, ranges));
    }

    double ihh1=0;
    double ihh0=0;

    int n_c0 = hm->hapData->get_n_c0(locus);
    int n_c1 = hm->hapData->get_n_c1(locus);

    double curr_ehh0_before_norm = 0;
    double curr_ehh1_before_norm = 0;

    double prev_ehh0_before_norm = 0;
    double prev_ehh1_before_norm = 0;

    const double& normalizer_0 = twice_num_pair_or_square(n_c0, p.ALT);
    const double& normalizer_1 = twice_num_pair_or_square(n_c1, p.ALT);

    int numSnps = hm->hapData->nloci;
    int numHaps = hm->hapData->nhaps;

    // PHASE 1: INITIALIZATION
    int* group_count = new int[numHaps];
    int* group_id = new int[numHaps];
    bool* isDerived = new bool [numHaps];
    // bool* isAncestral = new bool [numHaps];

    for(int i = 0; i<numHaps; i++){ //will be vectorized with compile time flags
        group_count[i] = 0;
        group_id[i] = 0;
        isDerived[i] = false;
        // isAncestral[i] = false;        //assert(hm->hapData->hapEntries[locus].flipped == false);
    }
    int totgc = 0; // total group count upto this point


    // /*** skip due to MAF ***/
    // if( hm->hapData->get_maf(locus) < p.MAF){
    //     {std::lock_guard<std::mutex> lock(mutex_log);
    //     (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
    //             << " (pos " << hm->mapData->mapEntries[locus].physicalPos << ") has MAF < " << p.MAF
    //             << ". Skipping calculation at this locus.\n";
    //     }//unlock
    //     return skipLocusPair();
    // }

    /*** skip due to MAC ***/
    if(n_c0 == 1 || n_c1 == 1){
        {std::lock_guard<std::mutex> lock(mutex_log);
        (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
                << " (pos " << hm->mapData->mapEntries[locus].physicalPos << ") has minor allele count = 1. Skipping calculation at this locus.\n";
        }//unlock
        return skipLocusPair();   
    }

    // PHASE 1a: INIT CORE LOCUS
    /* if monomorphic core was allowed, i.e. missing*/
    // if(n_c1==0){    // all 0s
    //     group_count[0] = numHaps;
    //     totgc+=1;
    //     curr_ehh0_before_norm = normalizer_0;
    // }else if (n_c1==numHaps){ // all 1s
    //     group_count[0] = numHaps;
    //     totgc+=1;
    //     MyBitset* vb = hm->hapData->hapEntries[locus].hapbitset;
    //     ACTION_ON_ALL_SET_BITS(vb, {
    //         isDerived[set_bit_pos] = true;
    //     });
    //     curr_ehh1_before_norm = normalizer_1;
    // }


    if(n_c1 == 0 || n_c0 == 0){ // monomorphic site
        if(downstream){ // assuming downstream is executed first
            std::lock_guard<std::mutex> lock(mutex_log); // lock the mutex
            (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
                << " (pos " << hm->mapData->mapEntries[locus].physicalPos << ") is monomorphic. Skipping calculation at this locus.\n"; // unlock happens automatically
        }
        return skipLocusPair();
    }else{  //so both n_c1 and n_c0 is non-0
        group_count[1] = n_c1;
        group_count[0] = n_c0;
        ACTION_ON_ALL_SET_BITS(hm->hapData->hapEntries[locus].hapbitset, {
            isDerived[set_bit_pos] = true;
            group_id[set_bit_pos] = 1;
        });
        totgc+=2;
        curr_ehh0_before_norm = normalizer_0;
        curr_ehh1_before_norm = normalizer_1;
    }

    // PHASE 2: IHS LOOP
    
    // if(downstream){
    //     if(normalizer_1!=0){
    //         ihh1 += (curr_ehh1_before_norm + prev_ehh1_before_norm) * 0.5 / normalizer_1;
    //     }
    //     if(normalizer_0!=0){
    //         ihh0 += (curr_ehh0_before_norm + prev_ehh0_before_norm) * 0.5 /  normalizer_0;
    //     }
    // }

    prev_ehh1_before_norm = curr_ehh1_before_norm;
    prev_ehh0_before_norm = curr_ehh0_before_norm;

    if(OUT_EHHS){
        if(downstream){
            ehhs_out.open(ehh_filename, std::ios::out);
            ehhs_out << std::fixed << "pdist\tgdist\tderEHH\tancEHH\tEHH" << "\n";
        }else{
            ehhs_out.open(ehh_filename, std::ios::out | std::ios::app);
        }
    }

    int i = locus;  
    int prev_index = locus; // for multi-MAF
    int skipped_due_to_multimaf = 0; // for multi-MAF
    while(true){ // Upstream: for ( int i = locus+1; i<all_positions.size(); i++ )

        bool edgeBreak = false;
        edgeBreak = (downstream)? (i-1 < 0) : (i+1 >= numSnps);
        //edgeBreak = nextLocOutOfBounds(i, downstream);
        if(edgeBreak) {
            {std::lock_guard<std::mutex> lock(mutex_log);
            (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << p.EHH_CUTOFF
                    << ". position: "<< hm->mapData->mapEntries[locus].physicalPos << " id: " << hm->mapData->mapEntries[locus].locusName << endl;
            if (!p.TRUNC){
                (*flog) << "Skipping calculation.";
            }
            (*flog) << endl;
            }//unlock
            if (!p.TRUNC){
                return skipLocusPair();
            }
            break;
        }
        
        i = (downstream)? i-1 : i+1; // update i 
        if(p.MULTI_MAF){
            if(hm->hapData->get_maf(i) < p.MAF){
                skipped_due_to_multimaf++;
                continue;
            }
        }

        assert(i >= 0 and i < numSnps);

        double distance =  geneticDistance(i, prev_index, downstream);
         if(p.CALC_NSL && !p.CALC_IHS){
            //distance = snpDistance(i, prev_index, downstream);
            distance = 1;
        }
        double scale = double(p.SCALE_PARAMETER) / double(physicalDistance(i, prev_index, downstream) );
        if (scale > 1) scale = 1;
        distance *= scale;
        prev_index = i; // save previous index for genetic distance calculation

        if(p.CALC_IHS && !p.CALC_NSL){
            if(curr_ehh1_before_norm*1.0/normalizer_1 <= p.EHH_CUTOFF and curr_ehh0_before_norm*1.0/normalizer_0  <= p.EHH_CUTOFF){   // or cutoff, change for benchmarking against hapbin
                //DBG("Break reason for locus "<<locus<<":: EHH_CUTOFF. "<<curr_ehh1_before_norm*1.0/normalizer_1<< "<=" << p.EHH_CUTOFF <<endl);
                break;
            }
        }

        if(hm->hapData->get_n_c0(i) == 0 or hm->hapData->get_n_c1(i) == 0 ){ 
            // monomorphic check, it is possible to skip compute to speed it up
            assert(normalizer_1!=0);
            ihh1 += (prev_ehh1_before_norm + curr_ehh1_before_norm) * distance * 0.5 / normalizer_1;
            
            assert(normalizer_0!=0);
            ihh0 += (prev_ehh0_before_norm + curr_ehh0_before_norm) * distance * 0.5 / normalizer_0;

            if(OUT_EHHS){
                string neg = downstream? "-" : ""; // if downstream, then negative distance
                double ehh1 = curr_ehh1_before_norm / normalizer_1;
                double ehh0 = curr_ehh0_before_norm / normalizer_0;
                double ehh = (curr_ehh1_before_norm + curr_ehh0_before_norm) / twice_num_pair_or_square(n_c1+n_c0, p.ALT);
                //if in range
                ehhs_out << std::fixed << neg << physicalDistance(i, locus, downstream)  << "\t"
                <<  neg <<  geneticDistance(i, locus, downstream) << "\t"
                << ehh1 << "\t"
                << ehh0 << "\t"
                << ehh << "";
                ehhs_out << endl;
            }
            continue;
        }

        // PHASE 2: GET MAP
        std::unique_ptr<std::unordered_map<int, std::vector<int>>> mp(new std::unordered_map<int, std::vector<int>>());
        unordered_map<int, vector<int> >& m = (* mp);

        MyBitset* v2p = hm->hapData->hapEntries[i].hapbitset;
        ACTION_ON_ALL_SET_BITS(v2p, {
            int old_group_id = group_id[set_bit_pos];
            m[old_group_id].push_back(set_bit_pos);
        });
        

        // PHASE 3: UPDATE EHH FROM SPLIT

        for (const auto &ele : m) {
            int old_group_id = ele.first;
            int newgroup_size = ele.second.size() ;
                            
            if(group_count[old_group_id] == newgroup_size || newgroup_size == 0){
                continue;
            }

            for(const int &v: ele.second){
                group_id[v] = totgc;
            }
            
            double del_update = -twice_num_pair_or_square(group_count[old_group_id], p.ALT) + twice_num_pair_or_square(newgroup_size, p.ALT) + twice_num_pair_or_square(group_count[old_group_id] - newgroup_size, p.ALT);
            
            group_count[old_group_id] -= newgroup_size;
            group_count[totgc] += newgroup_size;
            
            totgc+=1;
            
            // uncomment for flipped
            //bool isDerivedGroup =  (!hm->hapData->hapEntries[locus].flipped && isDerived[ele.second[0]]) || (hm->hapData->hapEntries[locus].flipped && !isAncestral[ele.second[0]]); // just check first element to know if it is derived. 
            
            bool isDerivedGroup =  isDerived[ele.second[0]];
            if(isDerivedGroup) // if the core locus for this chr has 1 (derived), then update ehh1, otherwise ehh0
            {
                curr_ehh1_before_norm += del_update;
            }else{
                curr_ehh0_before_norm += del_update;
            }
        }

        m.clear(); // CLEAR THE MAP //unordered_map<int, vector<int> >().swap(m);
        
        if(prev_ehh1_before_norm*1.0/normalizer_1 > p.EHH_CUTOFF){ // bug fix for low MAC
            assert(normalizer_1!=0);
            ihh1 += (prev_ehh1_before_norm + curr_ehh1_before_norm) * distance * 0.5 / normalizer_1;
            prev_ehh1_before_norm = curr_ehh1_before_norm;
        }

        if(prev_ehh0_before_norm*1.0/normalizer_0 > p.EHH_CUTOFF){ // bug fix for low MAC
            assert(normalizer_0!=0);
            ihh0 += (prev_ehh0_before_norm + curr_ehh0_before_norm) * distance * 0.5 / normalizer_0;
            prev_ehh0_before_norm = curr_ehh0_before_norm;
        }

        if(OUT_EHHS){
            string neg = downstream? "-" : ""; // if downstream, then negative distance
            double ehh1 = curr_ehh1_before_norm / normalizer_1;
            double ehh0 = curr_ehh0_before_norm / normalizer_0;
            double ehh = (curr_ehh1_before_norm + curr_ehh0_before_norm) / twice_num_pair_or_square(n_c1+n_c0, p.ALT);
            //if in range
            ehhs_out << std::fixed << neg << physicalDistance(i, locus, downstream)  << "\t"
            <<  neg <<  geneticDistance(i, locus, downstream) << "\t"
            << ehh1 << "\t"
            << ehh0 << "\t"
            << ehh << "";
            ehhs_out << endl;
        }


        if(totgc == numHaps) {
            DBG("Break reason for locus "<<locus<<":: ALL_UNIQUE."<<endl);
            break;
        }
        if(!p.CALC_NSL && p.CALC_IHS && physicalDistance(i,locus, downstream) >= max_extend) {
            DBG("Break reason for locus "<<locus<<":: MAX_EXTEND."<<endl);
            break;
        }
        if(p.CALC_NSL && !p.CALC_IHS && abs(i-locus)-skipped_due_to_multimaf >= max_extend) {
            DBG("Break reason for locus "<<locus<<":: MAX_EXTEND_NSL."<<endl);
            break; 
        }
    }




    delete[] group_count;
    delete[] group_id;
    delete[] isDerived;
    //delete[] isAncestral;
    return make_pair(ihh1, ihh0);
}



void IHS::main() {
    if(p.CALC_NSL && !p.CALC_IHS){
        this->max_extend = ( p.MAX_EXTEND_NSL <= 0 ) ? physicalDistance(0,hm->hapData->nloci-1,true) : p.MAX_EXTEND_NSL;
        init_global_fout("nsl");
    }else{
        this->max_extend = ( p.MAX_EXTEND <= 0 ) ? physicalDistance(0,hm->hapData->nloci-1,true) : p.MAX_EXTEND;
        init_global_fout("ihs");
    }

    DBG("DEBUG::: Total number of loci: "<<hm->hapData->nloci<<endl);

    // if(p.MISSING_ALLOWED){ // do imputation
    //     if(hm->p.MISSING_MODE == "NO_IMPUTE"){
    //         for(int locus = 0; locus <  hm->mapData->nloci; locus++) {  // TODO: make multithreaded
    //             pair<double, double> ihh1_ihh0 = infer_missing(locus);  
    //         }
    //         hm->hapData->assignVerdict();
    //     }
    // }

    std::thread progressBarThread(displayProgressBar, std::ref(hm->currentProcessed), hm->hapData->nloci); // Launch the progress bar in a separate thread
    ThreadPool pool(p.numThreads);

    if(!p.WRITE_DETAILED_IHS){
        string ihsornsl = p.CALC_NSL ? "nsl" : "ihs";
        if(hm->hapData->unphased){
            *fout << "id" << "\t" <<   "pos" << "\t"
                            << "freq" << "\t"
                            << "ihh1" << "\t" << "ihh0" << "\t"  << ihsornsl  << endl;
        }else{
            *fout << "id" << "\t" <<   "pos" << "\t"
                            << "freq" << "\t"
                            << "ihh1" << "\t" << "ihh0" << "\t"  << ihsornsl  << endl;
        }
        
        std::vector< std::future<pair<double, double> > > results;
        for(int i = 0; i <  hm->mapData->nloci; ++i) {

            //avoid doing this here,as Tasks may not be guaranteed to run correctly or complete in order.
                        // if(hm->mapData->mapEntries[i].locId == -1) { //p.MULTI_MAF
            //     continue;
            // }
            // if(hm->hapData->get_maf(i) < p.MAF) { // if core locus has MAF < p.MAF, skip it, useful in --keep-low-freq
            //     continue;
            // }

            results.emplace_back(
                pool.enqueue([i,this] {
                    return this->calc_ihh1(i);
                })
            );
            
        }

        int locus = 0; 
        for(auto && result: results){ // this is a blocking call
            pair<double, double> ihh1_ihh0 = result.get(); 
            double ihh1 = ihh1_ihh0.first;
            double ihh0 = ihh1_ihh0.second;

            if(!skipLocus(ihh1_ihh0)){
                if(hm->hapData->unphased){
                    const double& iHS2 = ihh1;
                    const double& iHS0 = ihh0;
                    double ihs = iHS2 > iHS0 ? iHS2 : 0-iHS0;
                    //std::fixed <<   std::setprecision(6) << 
                    *fout << hm->mapData->mapEntries[locus].locusName << "\t" <<   hm->mapData->mapEntries[locus].physicalPos << "\t"
                            << hm->hapData->calcFreq(locus) << "\t"
                            << iHS2 << "\t" << iHS0 <<"\t"<< ihs <<endl;
                            //<<  hm->mapData->mapEntries[locus].locId << "\t" 
                }else{  
                    //std::fixed <<   std::setprecision(6) <<  
                    *fout << hm->mapData->mapEntries[locus].locusName << "\t" <<   hm->mapData->mapEntries[locus].physicalPos << "\t"
                            << hm->hapData->calcFreq(locus) << "\t"
                            << ihh1 << "\t" << ihh0 <<"\t"<< log10(ihh1/ihh0) <<endl;
                            //hm->mapData->mapEntries[locus].locId << "\t"
                }
            }

            locus++;
            hm->currentProcessed++;
            assert(locus<=hm->mapData->nloci);
        }
        progressBarThread.join();
    }

    if(p.WRITE_DETAILED_IHS){
        string ihsornsl = p.CALC_NSL ? "nsl" : "ihs";
        if(hm->hapData->unphased){
            cerr<<"ERROR: Detailed iHS output not supported for unphased data. Please run without --write-detailed-ihs flag."<<endl;
            exit(1);
        }
        *fout << "id" << "\t" <<   "pos" << "\t"
                            << "freq" << "\t" << "ihh1" << "\t" << "ihh0" << "\t"  << ihsornsl << "\t"
                            << "ihh1_left" << "\t" << "ihh1_right" << "\t" << "ihh0_left" << "\t" << "ihh0_right" << endl;
        std::vector< std::future<IhhComponents> > results;
        for(int i = 0; i <  hm->mapData->nloci; ++i) {
            results.emplace_back(
                pool.enqueue([i,this] {
                    return this->calc_ihh1_details(i);
                })
            );            
        }

        int locus = 0; 
        //*fout << "#locusID\tphysPos\tfreq1\tihh1\tihh0\tunstd_iHS\tihh1_left\tihh1_right\tihh0_left\tihh0_right\n";
        for(auto && result: results){ // this is a blocking call
            IhhComponents ihh = result.get();

            // Only proceed if locus is valid
            if (!skipLocus(ihh.derived_right)) {
                // Total IHH values
                double ihh1 = ihh.derived_right + ihh.derived_left;
                double ihh0 = ihh.ancestral_right + ihh.ancestral_left;

                // Derived and ancestral IHH values    
                const auto& entry = hm->mapData->mapEntries[locus];
                double freq1 = hm->hapData->calcFreq(locus);

                // Compute unstandardized iHS value
                double unstandardized_iHS;
                if (hm->hapData->unphased) {
                    unstandardized_iHS = (ihh1 > ihh0) ? log10(ihh1 / ihh0) : -log10(ihh0 / ihh1); //ihs2==ihh1, ihs0==ihh0
                } else {
                    unstandardized_iHS = log10(ihh1 / ihh0);
                }

                *fout << entry.locusName << "\t"
                    << entry.physicalPos << "\t"
                    << freq1 << "\t"
                    << ihh1 << "\t"
                    << ihh0 << "\t"
                    << unstandardized_iHS << "\t"
                    << ihh.derived_left << "\t"
                    << ihh.derived_right << "\t"
                    << ihh.ancestral_left << "\t"
                    << ihh.ancestral_right << "\n";
            }
            locus++;
            hm->currentProcessed++;
            assert(locus<=hm->mapData->nloci);
        }
        progressBarThread.join();
    }

    DBG("DEBUG::: IHS done."<<endl);

    //output file closing handled by selscan-stats.h
}

    
/**
 * @brief Calculate the IHH statistics for a single locus (as part of iHS routine).
 *        Computes upstream and downstream integration of EHH decay for derived (1) and ancestral (0) alleles.
 * @param locus The locus index in the haplotype matrix.
 * @return Pair of iHS-like values: (derived_score, ancestral_score)
 */
std::pair<double, double> IHS::calc_ihh1(int locus) {

    if(hm->hapData->get_maf(locus) < p.MAF) { // if core locus has MAF < p.MAF, skip it, useful in --keep-low-freq
        {
                string reason = "MAF < " + std::to_string(p.MAF);
                std::lock_guard<std::mutex> lock(mutex_log);
                (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
                << " (pos " << hm->mapData->mapEntries[locus].physicalPos << ") has " 
                << reason << ". Skipping this locus." << std::endl;
            }
        return skipLocusPair();
    }

    // === Unphased mode ===
    if (hm->hapData->unphased) {
        // Compute downstream unphased iHH (derived and ancestral) and cumulative iHH (for normalization)
        OutputUnphasedIHH downstream = calc_ehh_unidirection_unphased(locus, true);  // true = downstream
        if (skipLocus(downstream.iHH0)) return skipLocusPair();

        // Compute upstream unphased iHH (derived and ancestral)
        OutputUnphasedIHH upstream = calc_ehh_unidirection_unphased(locus, false);  // false = upstream
        if (skipLocus(upstream.iHH0)) return skipLocusPair();

        // Combine totals
        double total_iHH2 = upstream.iHH2 + downstream.iHH2;
        double total_ciHH2 = upstream.ciHH2 + downstream.ciHH2;

        double total_iHH0 = upstream.iHH0 + downstream.iHH0;
        double total_ciHH0 = upstream.ciHH0 + downstream.ciHH0;

        if(total_iHH2 == 0 || total_ciHH2 == 0 || total_iHH0 == 0 || total_ciHH0 == 0){ 
            std::string reason;
            if (total_iHH2 == 0) reason = "iHH2 = 0";
            else if (total_ciHH2 == 0) reason = "ciHH2 = 0";
            else if (total_iHH0 == 0) reason = "iHH0 = 0";
            else if (total_ciHH0 == 0) reason = "ciHH0 = 0";

            {
                std::lock_guard<std::mutex> lock(mutex_log);
                (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
                << " (pos " << hm->mapData->mapEntries[locus].physicalPos << ") has " 
                << reason << ". Skipping this locus." << std::endl;
            }
            return skipLocusPair();
        }

        double ihs2 = log10(total_iHH2 / total_ciHH2);
        double ihs0 = log10(total_iHH0 / total_ciHH0);

        return { ihs2, ihs0 };
    }

    // === Phased mode ===

    pair<double, double> ihh1_ihh0_downstream = calc_ehh_unidirection(locus, true);
    if (skipLocus(ihh1_ihh0_downstream)) return skipLocusPair();

    pair<double, double> ihh1_ihh0_upstream = calc_ehh_unidirection(locus, false);
    if (skipLocus(ihh1_ihh0_upstream)) return skipLocusPair();

    // Sum values
    double ihh1 = ihh1_ihh0_upstream.first + ihh1_ihh0_downstream.first;
    double ihh0 = ihh1_ihh0_upstream.second + ihh1_ihh0_downstream.second;

    if (ihh1 == 0 || ihh0 == 0) {
        //should not happen, but just in case
        {
            string reason = "zero ihh";
                std::lock_guard<std::mutex> lock(mutex_log);
                (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
                << " (pos " << hm->mapData->mapEntries[locus].physicalPos << ") has " 
                << reason << ". Skipping this locus." << std::endl;
            }
        //DBG("WARNING: ihh is zero for locus " << locus << ". Likely due to early EHH decay or no derived homozygosity.");
        return skipLocusPair();
    }

    return { ihh1, ihh0 };
}


/**
 * @brief Calculate the IHH statistics for a single locus (as part of iHS routine).
 *        Computes right (upstream) and left (downstream) integration of EHH decay
 *        for derived (1) and ancestral (0) alleles.
 * @param locus The locus index in the haplotype matrix.
 * @return IhhComponents with directional values for derived and ancestral alleles.
 */
IhhComponents IHS::calc_ihh1_details(int locus) {
    if(hm->hapData->get_maf(locus) < p.MAF) { // if core locus has MAF < p.MAF, skip it, useful in --keep-low-freq
    {
        string reason = "MAF < " + std::to_string(p.MAF);
        std::lock_guard<std::mutex> lock(mutex_log);
        (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
        << " (pos " << hm->mapData->mapEntries[locus].physicalPos << ") has " 
        << reason << ". Skipping this locus." << std::endl;
        }
        return {SKIP_LOCUS_VALUE, SKIP_LOCUS_VALUE, SKIP_LOCUS_VALUE, SKIP_LOCUS_VALUE};
    }

    if (hm->hapData->unphased) {
        OutputUnphasedIHH left = calc_ehh_unidirection_unphased(locus, true);     // left (downstream)
        if(left.iHH0 == SKIP_LOCUS_VALUE){
            return {SKIP_LOCUS_VALUE, SKIP_LOCUS_VALUE, SKIP_LOCUS_VALUE, SKIP_LOCUS_VALUE};
        }

        OutputUnphasedIHH right = calc_ehh_unidirection_unphased(locus, false); // right (upstream)
        if(right.iHH0 == SKIP_LOCUS_VALUE){
            return {SKIP_LOCUS_VALUE, SKIP_LOCUS_VALUE, SKIP_LOCUS_VALUE, SKIP_LOCUS_VALUE};
        }

        return { right.iHH2, left.iHH2, right.iHH0, left.iHH0 };
    }else{
        // Phased data
        double ihh1_right = 0, ihh1_left = 0;
        double ihh0_right = 0, ihh0_left = 0;

        auto ihh1_ihh0_left = calc_ehh_unidirection(locus, true);   // left (downstream)
        if (skipLocus(ihh1_ihh0_left)) return {SKIP_LOCUS_VALUE, SKIP_LOCUS_VALUE, SKIP_LOCUS_VALUE, SKIP_LOCUS_VALUE};

        auto ihh1_ihh0_right = calc_ehh_unidirection(locus, false); // right (upstream)
        if (skipLocus(ihh1_ihh0_right)) return {SKIP_LOCUS_VALUE, SKIP_LOCUS_VALUE, SKIP_LOCUS_VALUE, SKIP_LOCUS_VALUE};

     
        ihh1_right = ihh1_ihh0_right.first;
        ihh0_right = ihh1_ihh0_right.second;
        ihh1_left = ihh1_ihh0_left.first;
        ihh0_left = ihh1_ihh0_left.second;

        return { ihh1_right, ihh1_left, ihh0_right, ihh0_left };
    }
}