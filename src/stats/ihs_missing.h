#include "ihs.h"

//======= MISSING RELATED FUNCTIONS =======//
pair<double, double> IHS::infer_missing(int locus){
    if(hm->p.UNPHASED){
        double a, b;
        pair<double, double>  downstream = calc_ehh_unidirection_unphased_missing(locus, true, a, b); // downstream
        pair<double, double> upstream = calc_ehh_unidirection_unphased_missing(locus, false,a, b); // upstream
    }else{
        pair<double, double>  downstream = calc_ehh_unidirection_missing(locus, true); // downstream
        pair<double, double> upstream = calc_ehh_unidirection_missing(locus, false); // upstream
    }
    
    return make_pair(1, 1);
}

/**
 * Calculate EHH in only one direction until cutoff is hit - upstream or downstream
*/
pair<double, double> IHS::calc_ehh_unidirection_missing(int locus, bool downstream){

    double ihh1=0;
    double ihh0=0;

    int n_c0 = hm->hapData->get_n_c0(locus);
    int n_c1 = hm->hapData->get_n_c1(locus);
    int n_c_missing = hm->hapData->get_n_c_missing(locus);

    double curr_ehh0_before_norm = 0;
    double curr_ehh1_before_norm = 0;

    double prev_ehh0_before_norm = 0;
    double prev_ehh1_before_norm = 0;

    double normalizer_0 = twice_num_pair_or_square(n_c0, p.ALT);
    double normalizer_1 = twice_num_pair_or_square(n_c1, p.ALT);

    int numSnps = hm->hapData->nloci;
    int numHaps = hm->hapData->nhaps;
    int numHaps_core = hm->hapData->nhaps - n_c_missing;


    // PHASE 1: INITIALIZATION
    int* group_count = new int[numHaps];
    int* group_id = new int[numHaps];
    bool* isDerived = new bool [numHaps];
    bool* isMissing = new bool [numHaps];
    bool* skipHap = new bool [numHaps];
    

    for(int i = 0; i<numHaps; i++){ //will be vectorized with compile time flags
        group_count[i] = 0;
        group_id[i] = 0;
        isDerived[i] = false;
        isMissing[i] = false;
        skipHap[i] = false;

        // isAncestral[i] = false;        //assert(hm->hapData->hapEntries[locus].flipped == false);
    }
    int totgc = 0; // total group count upto this point

    // PHASE 1a: INIT CORE LOCUS
    if(n_c1==0){    // all 0s
        group_count[0] = numHaps;
        totgc+=1;
        curr_ehh0_before_norm = normalizer_0;
//        hm->mapData->mapEntries[locus].skipLocus = true;
        //return make_pair(-9999,-9999);

        MyBitset* vmissing = hm->hapData->hapEntries[locus].xorbitset;
        ACTION_ON_ALL_SET_BITS(vmissing, {
            isMissing[set_bit_pos] = true;
            skipHap[set_bit_pos] = true;
        });

    }else if (n_c0==0){ // all 1s
        group_count[0] = numHaps;
        totgc+=1;
        
        MyBitset* vb = hm->hapData->hapEntries[locus].hapbitset;
        ACTION_ON_ALL_SET_BITS(vb, {
            isDerived[set_bit_pos] = true;
        });

        MyBitset* vmissing = hm->hapData->hapEntries[locus].xorbitset;
        ACTION_ON_ALL_SET_BITS(vmissing, {
            isMissing[set_bit_pos] = true;
            skipHap[set_bit_pos] = true;
        });
        
        curr_ehh1_before_norm = normalizer_1;
        //hm->mapData->mapEntries[locus].skipLocus = true;
    }else{  //so both n_c1 and n_c0 is non-0
        group_count[1] = n_c1;
        group_count[0] = n_c0;
        // bool assign_one = false;
        // if(n_c1> n_c0){
        //     assign_one = true;
        // }

        
        ACTION_ON_ALL_SET_BITS(hm->hapData->hapEntries[locus].hapbitset, {
            isDerived[set_bit_pos] = true;
            group_id[set_bit_pos] = 1;
        });
        MyBitset* vmissing = hm->hapData->hapEntries[locus].xorbitset;
        ACTION_ON_ALL_SET_BITS(vmissing, {
            //if(assign_one){
                // isDerived[set_bit_pos] = true;
                // group_id[set_bit_pos] = 1;
            //}
            skipHap[set_bit_pos] = true;
            isMissing[set_bit_pos] = true;
        });
        if(false){//if(assign_one){
            n_c1 = n_c1 + n_c_missing;
                normalizer_0 = twice_num_pair_or_square(n_c0, p.ALT);
                normalizer_1 = twice_num_pair_or_square(n_c1, p.ALT);
        }else{
                //n_c0 = n_c0 + n_c_missing;
                normalizer_0 = twice_num_pair_or_square(n_c0, p.ALT);
                normalizer_1 = twice_num_pair_or_square(n_c1, p.ALT);
        }
        

        // if(hm->mapData->mapEntries[locus].skipLocus){
        //     {std::lock_guard<std::mutex> lock(mutex_log);
        //     (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
        //             << " (number " << locus << ") is monomorphic. Skipping calculation at this locus.\n";
        //     }//unlock
        //     return make_pair(0,0);
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

        if(!p.CALC_NSL && p.CALC_IHS){
            // if(curr_ehh1_before_norm*1.0/normalizer_1 <= p.EHH_CUTOFF and curr_ehh0_before_norm*1.0/normalizer_0  <= p.EHH_CUTOFF){   // or cutoff, change for benchmarking against hapbin
            //     //std::cerr<<"Break reason for locus "<<locus<<":: EHH_CUTOFF."<<endl;
            //     break;
            // }
        }
    
        bool edgeBreak = false;
        edgeBreak = nextLocOutOfBounds(i, downstream);
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
            //break;
        }
        i= (downstream)? i-1 : i+1;
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
            std::cerr << "ERROR: physical position not in ascending order.\n"; 
            exit(2);
        }
        if(p.CALC_NSL && !p.CALC_IHS){
            distance = 1;
        }
        double scale = double(p.SCALE_PARAMETER) / double(physicalDistance(i, downstream) );
        if (scale > 1) scale = 1;
        distance *= scale;

        if(hm->hapData->get_n_c0(i) == 0 or hm->hapData->get_n_c1(i) == 0 ){ // monomorphic check, do not compute unnecessarily
            // {std::lock_guard<std::mutex> lock(mutex_log);
            // std::cerr<<"ERROR: Monomorphic site should not exist."<<endl;
            // std::cerr<< hm->hapData->get_n_c0(i) <<" n_c0 at locus "<< i <<endl; 
            // std::cerr<< hm->hapData->get_n_c1(i) <<" n_c1 at locus "<< i<< endl; 
            // }//unlock
            // throw 0;
            
            //if you wish to continue anyway
            if(normalizer_1!=0){    // case where not all are 0s
                ihh1 += (prev_ehh1_before_norm + curr_ehh1_before_norm) * distance * 0.5 / normalizer_1;
            }
            if(normalizer_0!=0){   // case where not all are 1s
                ihh0 += (prev_ehh0_before_norm + curr_ehh0_before_norm) * distance * 0.5 / normalizer_0;
            }
            continue;
        }

        // PHASE 2: GET MAP
        std::unique_ptr<std::unordered_map<int, std::vector<int>>> mp(new std::unordered_map<int, std::vector<int>>());
        unordered_map<int, vector<int> >& m = (* mp);
        
        
            MyBitset* v2p = hm->hapData->hapEntries[i].hapbitset;
            ACTION_ON_ALL_SET_BITS(v2p, {
                if(!skipHap[set_bit_pos]){
                    int old_group_id = group_id[set_bit_pos];
                    m[old_group_id].push_back(set_bit_pos);
                }
            });
        

        // missing
        if(p.MISSING_MODE == "NO_IMPUTE"){
            
                unordered_map<int, vector<int> > mmiss;
                MyBitset* vmissing = hm->hapData->hapEntries[i].xorbitset;
                ACTION_ON_ALL_SET_BITS(vmissing, {
                    if(!skipHap[set_bit_pos]){
                        int old_group_id = group_id[set_bit_pos];
                        mmiss[old_group_id].push_back(set_bit_pos);
                    }
                });
                ACTION_ON_ALL_SET_BITS(vmissing, {
                    if(!skipHap[set_bit_pos]){
                        int old_group_id = group_id[set_bit_pos];
                        int missing_count_at_i = mmiss[old_group_id].size();
                        int biggroup_size = group_count[old_group_id]  ;
                        int newgroup_size = m[old_group_id].size() ;  // ones assign to 1
                        int split_old_group_size = biggroup_size - newgroup_size - missing_count_at_i; // ones assign to 0, rest are missing
                        if(split_old_group_size < 0){
                        //    cout<<set_bit_pos<<" "<<i<<" "<<locus<<" "<<old_group_id<<" "<<biggroup_size<<" "<<newgroup_size<<" "<<missing_count_at_i<<endl;
                        }
                        assert(split_old_group_size >= 0);
                        int num_samples = newgroup_size  + split_old_group_size;
                        double p1 = newgroup_size/(1.0*num_samples);
                        double p0 = split_old_group_size/(1.0*num_samples);
                        int len = abs(locus - i);
                        if(num_samples > 1){
                            hm->hapData->insert_into_missing_matrix(set_bit_pos, i, len, num_samples, p0, p1);
                            //lock_guard<mutex> lock(hm->mapData->map_mutex);
                            //cout<< set_bit_pos<<" "<<i<<" "<<len<<" " << locus<<" NS"<<num_samples<<" "<<p0<<" "<<p1<<" " << split_old_group_size<<endl;
                            skipHap[set_bit_pos] = true;
                        }
                        
                        
                    }
                });
                // for (const auto &ele : mmiss) {
                //     m[ele.first] = ele.second;
                // }
                mmiss.clear();
            
        }
        

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
            if(normalizer_1!=0){
                ihh1 += (prev_ehh1_before_norm + curr_ehh1_before_norm) * distance * 0.5 / normalizer_1;
            }
            prev_ehh1_before_norm = curr_ehh1_before_norm;
        }

        if(curr_ehh0_before_norm*1.0/normalizer_0 > p.EHH_CUTOFF){
            if(normalizer_0!=0){
                ihh0 += (prev_ehh0_before_norm + curr_ehh0_before_norm) * distance * 0.5 / normalizer_0;
            }
            prev_ehh0_before_norm = curr_ehh0_before_norm;
        }

        if(totgc == numHaps - n_c_missing) {
            //std::cerr<<"Break reason for locus "<<locus<<":: ALL_UNIQUE."<<endl;
            break;
        }
        if(!p.CALC_NSL && p.CALC_IHS && physicalDistance_from_core(i,locus, downstream) >= max_extend) {
            //std::cerr<<"Break reason for locus "<<locus<<":: MAX_EXTEND."<<endl;
            break;
        }
        if(p.CALC_NSL && !p.CALC_IHS && abs(i-locus) >= max_extend) {
            //std::cerr<<"Break reason for locus "<<locus<<":: MAX_EXTEND_NSL."<<endl;
            break; 
        }
    }

    delete[] group_count;
    delete[] group_id;
    delete[] isDerived;
    //delete[] isAncestral;
    return make_pair(ihh1, ihh0);
}



pair<double, double> IHS::calc_ehh_unidirection_unphased_missing(int locus, bool downstream,  double& cihh2, double& cihh0){

    std::unique_ptr<std::unordered_map<int, std::vector<int>>> mp(new std::unordered_map<int, std::vector<int>>());
    unordered_map<int, vector<int> >& m = (* mp);

    int numSnps = hm->hapData->nloci;
    int numHaps = hm->hapData->nhaps;

    double prev_ehh_before_norm[3];
    double curr_ehh_before_norm[3];
    double  prev_cehh_before_norm[3]; 
    double  curr_cehh_before_norm[3];

    int n_c[3] = {0,0,0};
    int n_c_missing = 0;

    // int* group_count = new int[numHaps];
    // int* group_id =  new int[numHaps];
    // int* group_core =  new int[numHaps];
    // bool* is1 =  new bool[numHaps];
    // bool* is2 =  new bool[numHaps];
    // bool* isMissing = new bool[numHaps];
    // bool* skipHap = new bool[numHaps];

    int group_count[numHaps];
    int group_id[numHaps];
    int group_core[numHaps];
    bool is1[numHaps];
    bool is2[numHaps];
    bool isMissing[numHaps];
    bool skipHap[numHaps];

    double iHH[3] = {0,0,0};
    double ciHH[3] = {0,0,0};
    
    for(int i = 0; i<numHaps; i++){ //hopefully this will be vectorized with compile time flags
        group_count[i] = 0;
        is1[i] = false;
        is2[i] = false;
        isMissing[i] = false;
        skipHap[i] = false;
    }

    int totgc=0;

    n_c[1] = hm->hapData->get_n_c1(locus); //hapEntries[locus].hapbitset->num_1s;
    n_c[2] = hm->hapData->get_n_c2(locus); // hm->hapData->hapEntries[locus].xorbitset->num_1s;
    n_c_missing = hm->hapData->get_n_c_missing(locus);
    n_c[0] = numHaps - n_c[1] - n_c[2] - n_c_missing;

    //cerr<<"count:"<<n_c[0]<<" "<<n_c[1]<<" "<<n_c[2]<<" "<<n_c_missing<<endl;


    if(n_c[1] + n_c[2] + n_c[0] + n_c_missing!= numHaps){
        cerr<<"ERROR: n_c1 + n_c2 + n_c0 != numHaps"<<endl;
        exit(2);
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

    ACTION_ON_ALL_SET_BITS(hm->hapData->hapEntries[locus].missbitset, {
        isMissing[set_bit_pos] = true;
        group_id[set_bit_pos] = 0;
        skipHap[set_bit_pos] = true;
    });



    if(group_count[0] == numHaps - n_c_missing){ //monomorphic site
        totgc+=1;
    }else if(group_count[2] == 0 ){ //second clause is redundant
        if(group_count[0] + group_count[1] != numHaps - n_c_missing){
            cerr<<"ERROR: gc2==0 locus"<<locus<<endl;
            exit(2);
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
                << " (number " << locus << ") is monomorphic. Skipping calculation at this locus. "
                << "het: " <<  n_c[1] << " hom0: " << n_c[0] << " hom1: " << n_c[2] << ".\n";
        }
        //hm->mapData->mapEntries[locus].skipLocus = true;
        return skipLocusPair();
    }

    // double freqHetGT = n_c[1]*1.0/numHaps;
    // if (  freqHetGT > 1-p.MAF ) 
    // {
    //     {std::lock_guard<std::mutex> lock(mutex_log);
    //     (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
    //             << " (number " << locus << ") has too many hets. Skipping calculation at this locus. "
    //             << "het: " <<  n_c[1] << " hom0: " << n_c[0] << " hom1: " << n_c[2] << ".\n";
    //     }//unlock
    //     hm->mapData->mapEntries[locus].skipLocus = true;
    //     return make_pair(0,0);
    // }

    int i = locus;  // locus == core_locus
    while(true){ // Upstream: for ( int i = locus+1; i<all_positions.size(); i++ )
    
        
        // if(p.CALC_IHS && !p.CALC_NSL){
        //     if(curr_ehh_before_norm[2]*1.0/normalizer[2] <= p.EHH_CUTOFF and curr_ehh_before_norm[0]*1.0/normalizer[0]  <= p.EHH_CUTOFF){   // or cutoff, change for benchmarking against hapbin
        //         //std::cerr<<"Break reason for locus "<<locus<<":: EHH_CUTOFF."<<endl;
        //         break;
        //     }
        // }

        bool edgeBreak = false;
        edgeBreak = nextLocOutOfBounds(i, downstream); //(downstream)? (i-1 < 0) : (i+1 >= numSnps);
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
            std::cerr << "ERROR: physical position not in ascending order.\n"; 
            exit(2);
        }
        if(p.CALC_NSL && !p.CALC_IHS){
            distance = 1;
        }
        double scale = double(p.SCALE_PARAMETER) / double(physicalDistance(i, downstream) );
        if (scale > 1) scale = 1;
        distance *= scale;


        // if(i==22){
        //             cout<<"---"<<physicalDistance(i, downstream)<<endl;
        //             cout<<"---"<<geneticDistance(i, downstream)<<endl;
        //         }
        

        if (physicalDistance(i,downstream) > p.MAX_GAP)
        {
            {std::lock_guard<std::mutex> lock(mutex_log);
            (*flog) << "WARNING: Reached a gap of " << physicalDistance(i,downstream)
                    << "bp > " << p.MAX_GAP << "bp. Skipping calculation at position " << hm->mapData->mapEntries[locus].physicalPos << " id: " << hm->mapData->mapEntries[locus].locusName << "\n";
            }
            return skipLocusPair();
            //break;
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



        //if(p.MISSING_MODE == "NO_IMPUTE"){
        unordered_map<int, vector<int> > mmiss;
        unordered_map<int, vector<int> > m2;

        

        MyBitset* vmissing = hm->hapData->hapEntries[i].missbitset;
        ACTION_ON_ALL_SET_BITS(vmissing, {
            if(!skipHap[set_bit_pos]){
                int old_group_id = group_id[set_bit_pos];
                mmiss[old_group_id].push_back(set_bit_pos);
            }
        });


        //0-2 split
        ACTION_ON_ALL_SET_BITS(hm->hapData->hapEntries[i].xorbitset, {
            if(!skipHap[set_bit_pos]){
                int old_group_id = group_id[set_bit_pos];
                m2[old_group_id].push_back(set_bit_pos);
            }
        });
        

        //0-1 split
        ACTION_ON_ALL_SET_BITS(hm->hapData->hapEntries[i].hapbitset, {
            if(!skipHap[set_bit_pos]){
                int old_group_id = group_id[set_bit_pos];
                m[old_group_id].push_back(set_bit_pos);
            }
        });

        

        ACTION_ON_ALL_SET_BITS(vmissing, {
            if(!skipHap[set_bit_pos]){
                int old_group_id = group_id[set_bit_pos];
                int missing_count_at_i = mmiss[old_group_id].size();
                int biggroup_size = group_count[old_group_id]  ;

                int newgroup_size2 = m2[old_group_id].size() ;  // ones assign to 2
                int newgroup_size = m[old_group_id].size() ;  // ones assign to 1

                int split_old_group_size = biggroup_size - newgroup_size - newgroup_size2 - missing_count_at_i; // ones assign to 0, rest are missing

                //int split_old_group_size = biggroup_size - newgroup_size - missing_count_at_i; // ones assign to 0, rest are missing
                if(split_old_group_size < 0){
                    cout<<"DEBUG::: ncs" << hm->hapData->get_n_c0(i)<<" "<<hm->hapData->get_n_c1(i)<<" "<<hm->hapData->get_n_c2(i)<<" "<<hm->hapData->get_n_c_missing(i)<<endl;
                    cout<<"DEBUG::: "<<set_bit_pos<<" "<<i<<" "<<locus<<" "<<old_group_id<<" "<<biggroup_size<<" "<<newgroup_size<<" "<<newgroup_size2<<" "<<missing_count_at_i<<endl;
                }
                assert(split_old_group_size >= 0);
                int num_samples = newgroup_size  + newgroup_size2 + split_old_group_size; // non missing samples to take decision on
                double p0 = split_old_group_size/(1.0*num_samples);
                double p1 = newgroup_size/(1.0*num_samples);
                double p2 = newgroup_size2/(1.0*num_samples);
                int len = abs(locus - i);
                if(num_samples > 1){
                    hm->hapData->insert_into_missing_matrix(set_bit_pos, i, len, num_samples, p0, p1, p2);
                    //lock_guard<mutex> lock(hm->mapData->map_mutex);
                    //cout<< set_bit_pos<<" "<<i<<" "<<len<<" " << locus<<" NS"<<num_samples<<" "<<p0<<" "<<p1<<" " << split_old_group_size<<endl;
                    skipHap[set_bit_pos] = true;
                }
                
                
            }
        });

        

        for (const auto &ele : m2) {
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
        //updateEHH_from_split_unphased(m2, group_count, group_id, totgc, curr_ehh_before_norm, curr_cehh_before_norm, is1, is2, group_core);
        m2.clear();
        
        
        
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
        //updateEHH_from_split_unphased(m, group_count, group_id, totgc, curr_ehh_before_norm, curr_cehh_before_norm, is1, is2, group_core);
        m.clear();


        // missing
        if(p.MISSING_MODE == "NO_IMPUTE"){
            
            ACTION_ON_ALL_SET_BITS(vmissing, {
                if(!skipHap[set_bit_pos]){
                    int old_group_id = group_id[set_bit_pos];
                    int missing_count_at_i = mmiss[old_group_id].size();
                    int biggroup_size = group_count[old_group_id]  ;
                    int newgroup_size = m[old_group_id].size() ;  // ones assign to 1
                    int split_old_group_size = biggroup_size - newgroup_size - missing_count_at_i; // ones assign to 0, rest are missing
                    if(split_old_group_size < 0){
                    //    cout<<set_bit_pos<<" "<<i<<" "<<locus<<" "<<old_group_id<<" "<<biggroup_size<<" "<<newgroup_size<<" "<<missing_count_at_i<<endl;
                    }
                    assert(split_old_group_size >= 0);
                    int num_samples = newgroup_size  + split_old_group_size;
                    double p0 = newgroup_size/(1.0*num_samples);
                    double p1 = split_old_group_size/(1.0*num_samples);
                    int len = abs(locus - i);
                    if(num_samples > 1){
                        hm->hapData->insert_into_missing_matrix(set_bit_pos, i, len, num_samples, p0, p1);
                        //lock_guard<mutex> lock(hm->mapData->map_mutex);
                        //cout<< set_bit_pos<<" "<<i<<" "<<len<<" " << locus<<" NS"<<num_samples<<" "<<p0<<" "<<p1<<" " << split_old_group_size<<endl;
                        skipHap[set_bit_pos] = true;
                    }
                    
                    
                }
            });
            // for (const auto &ele : mmiss) {
            //     m[ele.first] = ele.second;
            // }
            mmiss.clear();
    
        }
        

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
            //std::cerr<<"Break reason for locus "<<locus<<":: ALL_UNIQUE."<<endl;
            break;
        }
        if ((!p.CALC_NSL && p.CALC_IHS) && physicalDistance_from_core(i, locus, downstream) >= max_extend) break;
        if ((p.CALC_NSL && !p.CALC_IHS) && abs(i-locus) >= max_extend) break;

        // if(downstream){
        //     cout<<locus<<":::l "<<i << " "<<curr_cehh_before_norm[0]/normalizer_not[0]<<" "<<curr_ehh_before_norm[0]/normalizer[0]<<" "<<ciHH[0]<<" "<<iHH[0]<<endl;

        // }else{
        //     cout<<locus<<":::r "<<i << " "<<curr_cehh_before_norm[0]/normalizer_not[0]<<" "<<curr_ehh_before_norm[0]/normalizer[0]<<" "<<ciHH[0]<<" "<<iHH[0]<<endl;
        // }
    }
    // delete[] group_count;
    // delete[] group_id;
    // delete[] is1;
    // delete[] is2;
    // delete[] group_core;
    // delete[] isMissing;
    // delete[] skipHap;
    

    cihh2 = ciHH[2];
    cihh0 = ciHH[0];
    return make_pair(iHH[2], iHH[0]);
}