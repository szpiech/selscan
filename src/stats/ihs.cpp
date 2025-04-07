#include "ihs.h"
#include <iomanip> 
#include <set>
#include <algorithm>
#include <cassert>
#include <pthread.h>


void IHS::updateEHH_from_split_unphased( unordered_map<int, vector<int> >& m, int* group_count, int* group_id, int& totgc, double* ehh_before_norm, double* cehh_before_norm, bool* is1, bool* is2, int* group_core){
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

pair<double, double> IHS::calc_ehh_unidirection_unphased(int locus, bool downstream,  double& cihh2, double& cihh0){

    std::unique_ptr<std::unordered_map<int, std::vector<int>>> mp(new std::unordered_map<int, std::vector<int>>());
    unordered_map<int, vector<int> >& m = (* mp);

    int numSnps = hm->hapData->nloci;
    int numHaps = hm->hapData->nhaps;

    double prev_ehh_before_norm[3];
    double curr_ehh_before_norm[3];
    double  prev_cehh_before_norm[3]; 
    double  curr_cehh_before_norm[3];

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
            HANDLE_ERROR("gc2==0 locus");
        }
        totgc+=2;
    }else{
        totgc+=3;
    }

    for (int i : {0, 2}) {
        //we dont need i = 1
        prev_ehh_before_norm[i] = curr_ehh_before_norm[i];
        prev_cehh_before_norm[i] =  curr_cehh_before_norm[i];
    }

    if(n_c[1] == numHaps || n_c[0] == numHaps || n_c[2] == numHaps){ 
        {
            std::lock_guard<std::mutex> lock(mutex_log);
            (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
            << " (pos " << hm->mapData->mapEntries[locus].physicalPos << ") is monomorphic. Skipping calculation at this locus. "
            << "het: " <<  n_c[1] << " hom0: " << n_c[0] << " hom1: " << n_c[2] << ".\n";
        }
        return skipLocusPair();
    }

    if(n_c[2] == 1){
        {std::lock_guard<std::mutex> lock(mutex_log);
        (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
                << " (pos " << hm->mapData->mapEntries[locus].physicalPos << ") has het count = 1. Skipping calculation at this locus.\n";
        }//unlock
        return skipLocusPair();
    }

    if(n_c[0] == 1){
        {std::lock_guard<std::mutex> lock(mutex_log);
        (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
                << " (pos " << hm->mapData->mapEntries[locus].physicalPos << ") has ancestral homozygotes count = 1. Skipping calculation at this locus.\n";
        }//unlock
        return skipLocusPair();
    }

    double freqHetGT = n_c[1]*1.0/numHaps;
    if (  freqHetGT > 1-p.MAF ) 
    {
        {
            std::lock_guard<std::mutex> lock(mutex_log);
            (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
                << " (pos " << hm->mapData->mapEntries[locus].physicalPos << ") has too many hets. Skipping calculation at this locus. "
                << "het: " <<  n_c[1] << " hom0: " << n_c[0] << " hom1: " << n_c[2] << ".\n";
        }//unlock
        return skipLocusPair();
    }

    int i = locus;  // locus == core_locus
    while(true){ // Upstream: for ( int i = locus+1; i<all_positions.size(); i++ )
    
        if(p.CALC_IHS && !p.CALC_NSL){
            if(curr_ehh_before_norm[2]*1.0/normalizer[2] <= p.EHH_CUTOFF and curr_ehh_before_norm[0]*1.0/normalizer[0]  <= p.EHH_CUTOFF){   // or cutoff, change for benchmarking against hapbin
                //DBG("Break reason for locus "<<locus<<":: EHH_CUTOFF.");
                break;
            }
        }

        bool edgeBreak = false;
        edgeBreak = nextLocOutOfBounds(i, downstream);
        if(edgeBreak) {
            {
                std::lock_guard<std::mutex> lock(mutex_log);
                (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << p.EHH_CUTOFF
                    << ". position: "<< hm->mapData->mapEntries[locus].physicalPos << " id: " << hm->mapData->mapEntries[locus].locusName << endl;
            }
            
            if (!p.TRUNC){
                {
                    std::lock_guard<std::mutex> lock(mutex_log);
                     (*flog) << "Skipping calculation at position " << hm->mapData->mapEntries[locus].physicalPos << " id: " << hm->mapData->mapEntries[locus].locusName << endl;
                }
                return skipLocusPair();
            }
            break;
        }
        i = (downstream) ? i-1 : i+1;
        
        
        // if(downstream){
        //     --i;
        //     if(i<0){
        //         break;
        //     }
        // }else{
        //     ++i;
        //     if(i>=numSnps){
        //         break;
        //     }
        // } 
        
        double distance =  geneticDistance(i, downstream);
        if (distance < 0) // this should not be needed as we already did integrity check previously
        {
            HANDLE_ERROR("physical position not in ascending order.");
        }
        if(p.CALC_NSL && !p.CALC_IHS){
            distance = 1;
        }
        double scale = double(p.SCALE_PARAMETER) / double(physicalDistance(i, downstream) );
        if (scale > 1) scale = 1;
        distance *= scale;


        if (physicalDistance(i,downstream) > p.MAX_GAP)
        {
            {
                std::lock_guard<std::mutex> lock(mutex_log);
                (*flog) << "WARNING: Reached a gap of " << physicalDistance(i,downstream)
                        << "bp > " << p.MAX_GAP << "bp. Skipping calculation at position " << hm->mapData->mapEntries[locus].physicalPos << " id: " << hm->mapData->mapEntries[locus].locusName << "\n";
            }
            return skipLocusPair();
        }

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

       if(hm->hapData->get_n_c0(i) == numHaps or hm->hapData->get_n_c1(i) == numHaps or hm->hapData->get_n_c2(i) == numHaps){ // monomorphic check, do not compute unnecessarily
            // {std::lock_guard<std::mutex> lock(mutex_log);
            // std::cerr<<"ERROR: Monomorphic site should not exist."<<endl;
            // std::cerr<< hm->hapData->get_n_c0(i) <<" n_c0 at locus "<< i <<endl; 
            // std::cerr<< hm->hapData->get_n_c1(i) <<" n_c1 at locus "<< i<< endl; 
            // }//unlock
            // throw 0;
            
            DBG("WARNING: Monomorphic site at locus "<<i);
            //if you wish to continue anyway
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
        
        ACTION_ON_ALL_SET_BITS(hm->hapData->hapEntries[i].xorbitset, {
            int old_group_id = group_id[set_bit_pos];
            m[old_group_id].push_back(set_bit_pos);
        });
        updateEHH_from_split_unphased(m, group_count, group_id, totgc, curr_ehh_before_norm, curr_cehh_before_norm, is1, is2, group_core);
        m.clear();

        ACTION_ON_ALL_SET_BITS(hm->hapData->hapEntries[i].hapbitset, {
            int old_group_id = group_id[set_bit_pos];
            m[old_group_id].push_back(set_bit_pos);
        });
        updateEHH_from_split_unphased(m, group_count, group_id, totgc, curr_ehh_before_norm, curr_cehh_before_norm, is1, is2, group_core);
        m.clear();

        // equivalent to calcHomozoygosity (without the normalization)
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

        if(double(curr_ehh_before_norm[0]*1.0/normalizer[0]) > p.EHH_CUTOFF || (p.CALC_NSL && !p.CALC_IHS) ){
            curr_ehh_before_norm[0] = ehh0_before_norm;
            if( normalizer[0]!=0){
                iHH[0] += ((curr_ehh_before_norm[0]/ normalizer[0])* 0.5 + (prev_ehh_before_norm[0]* 1.0/ normalizer[0]) * 0.5 ) *distance;
            }
            prev_ehh_before_norm[0] = curr_ehh_before_norm[0];
        }

        if(double(curr_ehh_before_norm[2]*1.0/normalizer[2]) > p.EHH_CUTOFF || (p.CALC_NSL && !p.CALC_IHS)){
            curr_ehh_before_norm[2] = ehh2_before_norm;
            if( normalizer[2]!=0){
                iHH[2] += ((curr_ehh_before_norm[2] / normalizer[2]) * 0.5 + (prev_ehh_before_norm[2] / normalizer[2]) * 0.5 ) *distance;
            }
            prev_ehh_before_norm[2] = curr_ehh_before_norm[2];
        }
            
        if(true){  //this is how it's in selscan, does not depend on cutoff
            curr_cehh_before_norm[0] = cehh0_before_norm;
            if(normalizer_not[0]!=0){
                ciHH[0] += ((curr_cehh_before_norm[0] / normalizer_not[0])*0.5  + (prev_cehh_before_norm[0]/ normalizer_not[0])*0.5)*distance ;
            }
            prev_cehh_before_norm[0] = curr_cehh_before_norm[0];

            curr_cehh_before_norm[2] = cehh2_before_norm;
            if(normalizer_not[2]!=0){
                ciHH[2] += ((curr_cehh_before_norm[2] / normalizer_not[2])*0.5  + (prev_cehh_before_norm[2] / normalizer_not[2])*0.5)*distance ;
            }
            prev_cehh_before_norm[2] = curr_cehh_before_norm[2];
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
            //DBG("Break reason for locus "<<locus<<":: ALL_UNIQUE.");
            break;
        }

        // if(i==16 && locus==172){
        //     cout<<locus<<":::l "<<i << " "<<curr_cehh_before_norm[0]/normalizer_not[0]<<" "<<curr_ehh_before_norm[0]/normalizer[0]<<" "<<ciHH[0]<<" "<<iHH[0]<<endl;
        //     cout<<hm->mapData->mapEntries[locus].physicalPos <<" "<<hm->mapData->mapEntries[i].physicalPos<<endl;
        // }

        if ((!p.CALC_NSL && p.CALC_IHS) && physicalDistance_from_core(i, locus, downstream) >= max_extend) {
            //DBG("Break reason for locus "<<locus<<":: MAX_EXTEND.");
            break;
        }
        if ((p.CALC_NSL && !p.CALC_IHS) && abs(i-locus) >= max_extend) {
            //DBG("Break reason for locus "<<locus<<":: MAX_EXTEND.");
            break;
        }

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

    cihh2 = ciHH[2];
    cihh0 = ciHH[0];
    return make_pair(iHH[2], iHH[0]);
}



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

    if(p.MULTI_MAF && hm->hapData->get_maf(locus) <= p.MAF){
        return skipLocusPair();
    }

    double ihh1=0;
    double ihh0=0;

    int prev_index = locus;
    int curr_index = -10000000;

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

    if(n_c0 == 1 || n_c1 == 1){
        {std::lock_guard<std::mutex> lock(mutex_log);
        (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
                << " (pos " << hm->mapData->mapEntries[locus].physicalPos << ") has minor allele count = 1. Skipping calculation at this locus.\n";
        }//unlock
        return skipLocusPair();
        
    }
    // if( hm->mapData->mapEntries[locus].locusName == "Locus0"){
    //     std::cerr<<"Locus0: "<<locus<<endl;
    //     std::cerr<<"n_c0: "<<n_c0<<endl;
    //     std::cerr<<"n_c1: "<<n_c1<<endl;
    //     std::cerr<<"normalizer_0: "<<normalizer_0<<endl;
    //     std::cerr<<"normalizer_1: "<<normalizer_1<<endl;
    // }
    // PHASE 1a: INIT CORE LOCUS
    if(n_c1==0){    // all 0s
        group_count[0] = numHaps;
        totgc+=1;
        curr_ehh0_before_norm = normalizer_0;
        //hm->mapData->mapEntries[locus].skipLocus = true;


        if(downstream){std::lock_guard<std::mutex> lock(mutex_log);
        (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
                << " (number " << locus << ") is monomorphic. Skipping calculation at this locus.\n";
        }//unlock
        
        return skipLocusPair();
    }else if (n_c1==numHaps){ // all 1s
        group_count[0] = numHaps;
        totgc+=1;
        MyBitset* vb = hm->hapData->hapEntries[locus].hapbitset;
        ACTION_ON_ALL_SET_BITS(vb, {
            isDerived[set_bit_pos] = true;
        });
        curr_ehh1_before_norm = normalizer_1;

        if(downstream){std::lock_guard<std::mutex> lock(mutex_log);
        (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
                << " (pos " << hm->mapData->mapEntries[locus].physicalPos << ") is monomorphic. Skipping calculation at this locus.\n";
        }//unlock
        return skipLocusPair();
        //hm->mapData->mapEntries[locus].skipLocus = true;
    }else{  //so both n_c1 and n_c0 is non-0
        group_count[1] = n_c1;
        group_count[0] = n_c0;
        ACTION_ON_ALL_SET_BITS(hm->hapData->hapEntries[locus].hapbitset, {
            isDerived[set_bit_pos] = true;
            group_id[set_bit_pos] = 1;
        });

        // if(hm->mapData->mapEntries[locus].skipLocus){
            // {std::lock_guard<std::mutex> lock(mutex_log);
            // (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
            //         << " (number " << locus << ") is monomorphic. Skipping calculation at this locus.\n";
            // }//unlock
            // return make_pair(0,0);
        // }
        
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

    int i = locus;  
    while(true){ // Upstream: for ( int i = locus+1; i<all_positions.size(); i++ )

        if(p. CALC_IHS && !p.CALC_NSL){
            if(curr_ehh1_before_norm*1.0/normalizer_1 <= p.EHH_CUTOFF and curr_ehh0_before_norm*1.0/normalizer_0  <= p.EHH_CUTOFF){   // or cutoff, change for benchmarking against hapbin
                //DBG("Break reason for locus "<<locus<<":: EHH_CUTOFF. "<<curr_ehh1_before_norm*1.0/normalizer_1<< "<=" << p.EHH_CUTOFF <<endl);
                break;
            }
        }
    
        bool edgeBreak = false;
        edgeBreak = (downstream)? (i-1 < 0) : (i+1 >= numSnps);
        //edgeBreak = nextLocOutOfBounds(i, downstream);
        if(edgeBreak) {
            {std::lock_guard<std::mutex> lock(mutex_log);
            (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << p.EHH_CUTOFF
                    << ". ";
            if (!p.TRUNC){
                //hm->mapData->mapEntries[locus].skipLocus = true;
                
                (*flog) << "Skipping calculation at position " << hm->mapData->mapEntries[locus].physicalPos << " id: " << hm->mapData->mapEntries[locus].locusName;
            }
            (*flog) << endl;
            }//unlock
            if (!p.TRUNC){
                return skipLocusPair();
            }
            break;
        }
        
        //update i 
        i = (downstream)? i-1 : i+1;
        assert(i>=0 and i<numSnps);

        if(p.MULTI_MAF){
            if(hm->hapData->get_maf(i) < p.MAF && p.SKIP){
                continue;
            }
        }
        // if(nextLocOutOfBounds(i, downstream)){
        //     break;
        // }

        // // if(edgeBreak) {
        // //     break;
        // // }
        // if(downstream){
        //     if (--i < 0) {
        //         //std::cerr<<"Break reason for locus "<<locus<<":: REACHED_LEFT_EDGE."<<endl;
        //         break;
        //     }
        // }else{
        //     if (++i >= numSnps) {
        //         //std::cerr<<"Break reason for locus "<<locus<<":: REACHED_RIGHT_EDGE."<<endl;
        //         break;
        //     }
        // }

        
        

        double distance =  geneticDistance(i, downstream);
        if (distance < 0) // this should not be needed as we already did integrity check previously
        {
            HANDLE_ERROR("physical position not in ascending order.");
        }
        double scale = double(p.SCALE_PARAMETER) / double(physicalDistance(i, downstream) );
        
        if(p.MULTI_MAF){
            curr_index = i;

            if(downstream){
                assert(curr_index < prev_index);
            }else{
                assert(curr_index > prev_index);
            }
            distance= abs(hm->mapData->mapEntries[prev_index].geneticPos - hm->mapData->mapEntries[curr_index].geneticPos);
            int phys_distance= abs(hm->mapData->mapEntries[prev_index].physicalPos - hm->mapData->mapEntries[curr_index].physicalPos);
            scale = double(p.SCALE_PARAMETER) / double(phys_distance);
            prev_index = curr_index;
        }

        if(p.CALC_NSL && !p.CALC_IHS){
            distance = 1;
        }
        if (scale > 1) scale = 1;
        distance *= scale;

        if(hm->hapData->get_n_c0(i) == 0 or hm->hapData->get_n_c1(i) == 0 ){ 
            
            // monomorphic check, can skip compute 
            //if you wish to update integral
            assert(normalizer_1!=0);
            ihh1 += (prev_ehh1_before_norm + curr_ehh1_before_norm) * distance * 0.5 / normalizer_1;
            
            assert(normalizer_0!=0);
            ihh0 += (prev_ehh0_before_norm + curr_ehh0_before_norm) * distance * 0.5 / normalizer_0;
            
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
        
        if(curr_ehh1_before_norm*1.0/normalizer_1 > p.EHH_CUTOFF){
            assert(normalizer_1!=0);
            ihh1 += (prev_ehh1_before_norm + curr_ehh1_before_norm) * distance * 0.5 / normalizer_1;
            prev_ehh1_before_norm = curr_ehh1_before_norm;
        }

        if(curr_ehh0_before_norm*1.0/normalizer_0 > p.EHH_CUTOFF){
            assert(normalizer_0!=0);
            ihh0 += (prev_ehh0_before_norm + curr_ehh0_before_norm) * distance * 0.5 / normalizer_0;
            prev_ehh0_before_norm = curr_ehh0_before_norm;
        }

        if(totgc == numHaps) {
            //DBG("Break reason for locus "<<locus<<":: ALL_UNIQUE."<<endl);
            break;
        }
        if(!p.CALC_NSL && p.CALC_IHS && physicalDistance_from_core(i,locus, downstream) >= max_extend) {
            //DBG("Break reason for locus "<<locus<<":: MAX_EXTEND."<<endl);
            break;
        }
        if(p.CALC_NSL && !p.CALC_IHS && abs(i-locus) >= max_extend) {
            //DBG("Break reason for locus "<<locus<<":: MAX_EXTEND_NSL."<<endl);
            break; 
        }
    }
    // if (ihh1 == 0) {
    //     DBG("ihh1 is zero for locus " << locus << ". Likely due to early EHH decay or no derived homozygosity.");
    // }
    // if (ihh0 == 0) {
    //     DBG("ihh0 is zero for locus " << locus << ". Likely due to early EHH decay or no ancestral homozygosity.");
    // }
    delete[] group_count;
    delete[] group_id;
    delete[] isDerived;
    //delete[] isAncestral;
    return make_pair(ihh1, ihh0);
}



void IHS::main() {
    if(p.CALC_NSL && !p.CALC_IHS){
        this->max_extend = ( p.MAX_EXTEND_NSL <= 0 ) ? physicalDistance_from_core(0,hm->hapData->nloci-1,true) : p.MAX_EXTEND_NSL;
        init_global_fout("nsl");
    }else{
        this->max_extend = ( p.MAX_EXTEND <= 0 ) ? physicalDistance_from_core(0,hm->hapData->nloci-1,true) : p.MAX_EXTEND;
        init_global_fout("ihs");
    }

    int total_calc_to_be_done = hm->hapData->nloci;
    DBG("DEBUG::: Total number of loci: "<<total_calc_to_be_done<<endl);

    
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
        std::vector< std::future<pair<double, double> > > results;
        for(int i = 0; i <  hm->mapData->nloci; ++i) {
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

    #ifdef DEBUG
        if(hm->p.MISSING_ALLOWED){
            hm->hapData->printMissingMatrix();
        }else{
            hm->hapData->print(5);
            //hm->hapData->print_by_ij();
        }
    #endif
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
    const int numSnps = hm->mapData->nloci;
    const int numHaps = hm->hapData->nhaps;

    // === Unphased mode ===
    if (hm->hapData->unphased) {
        double ihh2 = 0;
        double ihh0 = 0;
        double cihh2_upstream = 0;
        double cihh0_upstream  = 0;
        double cihh2_downstream = 0;
        double cihh0_downstream  = 0;

        pair<double, double>  ihh2_ihh0_downstream = calc_ehh_unidirection_unphased(locus, true,  std::ref(cihh2_downstream), std::ref(cihh0_downstream)); // downstream
        if (skipLocus(ihh2_ihh0_downstream)) return skipLocusPair();

        pair<double, double> ihh2_ihh0_upstream = calc_ehh_unidirection_unphased(locus, false, std::ref(cihh2_upstream), std::ref(cihh0_upstream)); // upstream
        if(skipLocus(ihh2_ihh0_upstream)) return skipLocusPair();

        ihh2 = ihh2_ihh0_upstream.first + ihh2_ihh0_downstream.first;
        ihh0 = ihh2_ihh0_upstream.second + ihh2_ihh0_downstream.second;

        double ihs2 = log10(ihh2/(cihh2_upstream+cihh2_downstream));
        double ihs0 = log10(ihh0/(cihh0_upstream+cihh0_downstream));

        if(ihh2 == 0) {
            {
                std::lock_guard<std::mutex> lock(mutex_log);
                (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
                << " (pos " <<  hm->mapData->mapEntries[locus].physicalPos << ") has ihh2 = 0. Skipping this locus. "<<endl;
            }
            return skipLocusPair();
        }

        if(ihh0 == 0){
            {
                std::lock_guard<std::mutex> lock(mutex_log);
                (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
                << " (pos " <<  hm->mapData->mapEntries[locus].physicalPos << ") has ihh0 = 0. Skipping this locus. "<<endl;
            }
            return skipLocusPair();
        }
        
        return { ihs2, ihs0 };
    }

    // === Phased mode ===
    double ihh1_up = 0, ihh1_down = 0;
    double ihh0_up = 0, ihh0_down = 0;

    pair<double, double> ihh1_ihh0_upstream = calc_ehh_unidirection(locus, false);
    if (skipLocus(ihh1_ihh0_upstream)) return skipLocusPair();

    pair<double, double> ihh1_ihh0_downstream = calc_ehh_unidirection(locus, true);
    if (skipLocus(ihh1_ihh0_downstream)) return skipLocusPair();

    // Sum values
    double ihh1 = ihh1_ihh0_upstream.first + ihh1_ihh0_downstream.first;
    double ihh0 = ihh1_ihh0_upstream.second + ihh1_ihh0_downstream.second;

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
    const int numSnps = hm->mapData->nloci;
    const int numHaps = hm->hapData->nhaps;

    if (hm->hapData->unphased) {
        double ihh2_right = 0, ihh2_left = 0;
        double ihh0_right = 0, ihh0_left = 0;

        // Calculate directional IHHs for unphased data
        calc_ehh_unidirection_unphased(locus, true, std::ref(ihh2_left), std::ref(ihh0_left));     // left (downstream)
        calc_ehh_unidirection_unphased(locus, false, std::ref(ihh2_right), std::ref(ihh0_right));  // right (upstream)

        return { ihh2_right, ihh2_left, ihh0_right, ihh0_left };
    }

    // Phased data
    double ihh1_right = 0, ihh1_left = 0;
    double ihh0_right = 0, ihh0_left = 0;

    auto ihh1_ihh0_right = calc_ehh_unidirection(locus, false); // right (upstream)
    if (skipLocus(ihh1_ihh0_right)) return {SKIP_LOCUS_VALUE, SKIP_LOCUS_VALUE, SKIP_LOCUS_VALUE, SKIP_LOCUS_VALUE};

    auto ihh1_ihh0_left = calc_ehh_unidirection(locus, true);   // left (downstream)
    if (skipLocus(ihh1_ihh0_left)) return {SKIP_LOCUS_VALUE, SKIP_LOCUS_VALUE, SKIP_LOCUS_VALUE, SKIP_LOCUS_VALUE};

    ihh1_right = ihh1_ihh0_right.first;
    ihh0_right = ihh1_ihh0_right.second;
    ihh1_left = ihh1_ihh0_left.first;
    ihh0_left = ihh1_ihh0_left.second;

    return { ihh1_right, ihh1_left, ihh0_right, ihh0_left };
}



// #include "ihs_missing.h"