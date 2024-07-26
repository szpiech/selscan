#include "ehh.h"

#define ACTION_ON_ALL_SET_BITS(hapbitset, ACTION)         \
    for (int k = 0; k < (hapbitset->nwords); k++) {             \
        uint64_t bitset = (hapbitset->bits)[k];                 \
        while (bitset != 0) {                        \
            uint64_t t = bitset & -bitset;           \
            int r = __builtin_ctzl(bitset);          \
            int set_bit_pos = (k * 64 + r);          \
            bitset ^= t;                             \
            ACTION;                                  \
        }                                            \
    }


/**
 * Calculate EHH in only one direction until cutoff is hit - upstream or downstream
*/
void EHH::calc_ehh_unidirection(int locus, bool downstream){

    int n_c0 = hm->hapData->get_n_c0(locus);
    int n_c1 = hm->hapData->get_n_c1(locus);

    double curr_ehh0_before_norm = 0;
    double curr_ehh1_before_norm = 0;
    double curr_ehh_before_norm = 0;

    const double& normalizer_0 = twice_num_pair_or_square(n_c0, p.ALT);
    const double& normalizer_1 = twice_num_pair_or_square(n_c1, p.ALT);
    const double& normalizer = twice_num_pair_or_square(n_c1+n_c0, p.ALT);


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

    // PHASE 1a: INIT CORE LOCUS
    if(n_c1==0){    // all 0s
        group_count[0] = numHaps;
        totgc+=1;
        curr_ehh0_before_norm = normalizer_0;
        hm->mapData->mapEntries[locus].skipLocus = true;
    }else if (n_c1==numHaps){ // all 1s
        group_count[0] = numHaps;
        totgc+=1;
        if(p.LOW_MEM){
            MyBitset* vb = hm->hapData->hapEntries[locus].hapbitset;
            ACTION_ON_ALL_SET_BITS(vb, {
                isDerived[set_bit_pos] = true;
            });
        }else{
            for (const int& set_bit_pos : hm->hapData->hapEntries[locus].positions){
                isDerived[set_bit_pos] = true;
            }
        }
        curr_ehh1_before_norm = normalizer_1;
        hm->mapData->mapEntries[locus].skipLocus = true;
    }else{  //so both n_c1 and n_c0 is non-0
        group_count[1] = n_c1;
        group_count[0] = n_c0;
        if(p.LOW_MEM){
            ACTION_ON_ALL_SET_BITS(hm->hapData->hapEntries[locus].hapbitset, {
                isDerived[set_bit_pos] = true;
                group_id[set_bit_pos] = 1;
            });
        }else{
            for (const int& set_bit_pos : hm->hapData->hapEntries[locus].positions){
                isDerived[set_bit_pos] = true;
                group_id[set_bit_pos] = 1;
            }
        }

        if(hm->mapData->mapEntries[locus].skipLocus){
            (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
                    << " (number " << locus << ") is monomorphic. \n";
        }
        
        totgc+=2;
        curr_ehh0_before_norm = normalizer_0;
        curr_ehh1_before_norm = normalizer_1;
    }
    curr_ehh_before_norm = normalizer_0 + normalizer_1;


    output[locus].ehh0 = curr_ehh0_before_norm*1.0/normalizer_0;
    output[locus].ehh1 = curr_ehh1_before_norm*1.0/normalizer_1;
    output[locus].ehh = curr_ehh_before_norm*1.0/normalizer;
    output[locus].gdist = 0;
    output[locus].pdist = 0;
    output[locus].print = true;
    

    int i = locus;  
    while(true){ // Upstream: for ( int i = locus+1; i<all_positions.size(); i++ )

        // if(curr_ehh1_before_norm*1.0/normalizer_1 <= p.EHH_CUTOFF and curr_ehh0_before_norm*1.0/normalizer_0  <= p.EHH_CUTOFF){   // or cutoff, change for benchmarking against hapbin
        //     //std::cerr<<"Break reason for locus "<<locus<<":: EHH_CUTOFF."<<endl;
        //     break;
        // }
    
        // bool edgeBreak = false;
        // edgeBreak = (downstream)? (i-1 < 0) : (i+1 >= numSnps);
        // if(edgeBreak) {
        //     (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << p.EHH_CUTOFF
        //             << ". ";
        //     (*flog) << endl;
        //     //break;
        // }
        
        if(downstream){
            if (--i < 0) {
                //std::cerr<<"Break reason for locus "<<locus<<":: REACHED_LEFT_EDGE."<<endl;
                break;
            }
        }else{
            if (++i >= numSnps) {
                //std::cerr<<"Break reason for locus "<<locus<<":: REACHED_RIGHT_EDGE."<<endl;
                break;
            }
        }

        double distance =  geneticDistance(i, downstream);
        if (distance < 0) // this should not be needed as we already did integrity check previously
        {
            std::cerr << "ERROR: physical position not in ascending order.\n"; 
            throw 0;
        }
        if(p.CALC_NSL){
            distance = 1;
        }
        double scale = double(p.SCALE_PARAMETER) / double(physicalDistance(i, downstream) );
        if (scale > 1) scale = 1;
        distance *= scale;

        if(hm->hapData->get_n_c0(i) == 0 or hm->hapData->get_n_c1(i) == 0 ){ // monomorphic check, do not compute unnecessarily
            // pthread_mutex_lock(&mutex_log);
            // std::cerr<<"ERROR: Monomorphic site should not exist."<<endl;
            // std::cerr<< hm->hapData->get_n_c0(i) <<" n_c0 at locus "<< i <<endl; 
            // std::cerr<< hm->hapData->get_n_c1(i) <<" n_c1 at locus "<< i<< endl; 
            // pthread_mutex_unlock(&mutex_log);
            // throw 0;
            
            //if you wish to continue anyway
            double current_derived_ehh = 0;
            double current_ancestral_ehh = 0;

            if(normalizer_1!=0){    // case where not all are 0s
                current_derived_ehh = curr_ehh1_before_norm*1.0/normalizer_1;
            }
            if(normalizer_0!=0){   // case where not all are 1s
                current_ancestral_ehh = curr_ehh0_before_norm*1.0/normalizer_0;
            }

            double current_ehh = curr_ehh_before_norm*1.0/normalizer;

            output[i].ehh0 = current_ancestral_ehh;
            output[i].ehh1 = current_derived_ehh;
            output[i].ehh = current_ehh;
            output[i].pdist = downstream? -physicalDistance_from_core(i,locus, downstream): physicalDistance_from_core(i,locus, downstream);
            output[i].gdist = downstream? -geneticDistance_from_core(i,locus, downstream): geneticDistance_from_core(i,locus, downstream);
            output[i].print = true;

            // if(downstream){
            //     (*fout) << std::fixed <<   -(hm->mapData->mapEntries[locus].physicalPos -  hm->mapData->mapEntries[i].physicalPos)  << "\t"
            //     <<  -(hm->mapData->mapEntries[locus].geneticPos -  hm->mapData->mapEntries[i].geneticPos)<< "\t"
            //     << current_derived_ehh << "\t"
            //     << current_ancestral_ehh << "\t";
            // }else{
            //     (*fout) << std::fixed <<   hm->mapData->mapEntries[i].physicalPos -  hm->mapData->mapEntries[locus].physicalPos  << "\t"
            //     <<  hm->mapData->mapEntries[i].geneticPos -  hm->mapData->mapEntries[locus].geneticPos<< "\t"
            //     << current_derived_ehh << "\t"
            //     << current_ancestral_ehh << "\t";
            // }
            // (*fout) << current_ehh << endl;


            continue;
        }

        // PHASE 2: GET MAP
        std::unique_ptr<std::unordered_map<int, std::vector<int>>> mp(new std::unordered_map<int, std::vector<int>>());
        unordered_map<int, vector<int> >& m = (* mp);

        if(p.LOW_MEM){
            MyBitset* v2p = downstream? hm->hapData->hapEntries[i+1].xorbitset : hm->hapData->hapEntries[i].xorbitset;
            // ensure that in boundary we don't do any xor calculation
            if( (hm->hapData->hapEntries[i].hapbitset->num_1s < v2p->num_1s && i!=numHaps-1) ||  hm->hapData->benchmark_flag != "XOR"){ 
                v2p = hm->hapData->hapEntries[i].hapbitset;
            }
            //v2p =hm->hapData->hapEntries[i].hapbitset; // uncomment to disable xor

            ACTION_ON_ALL_SET_BITS(v2p, {
                int old_group_id = group_id[set_bit_pos];
                m[old_group_id].push_back(set_bit_pos);
            });
        }else{
            vector<int> * v2p = downstream? &(hm->hapData->hapEntries[i+1].xors) : &(hm->hapData->hapEntries[i].xors);

             // ensure that in boundary we don't do any xor calculation
            if( (hm->hapData->hapEntries[i].positions.size() <(*v2p).size() && i!=numHaps-1) ||  hm->hapData->benchmark_flag != "XOR"){ 
                v2p = &hm->hapData->hapEntries[i].positions;
            }
            //v2p = &hm->hapData->hapEntries[i].positions; // uncomment to disable xor

            for (const int& set_bit_pos : (*v2p)){
                int old_group_id = group_id[set_bit_pos];
                m[old_group_id].push_back(set_bit_pos);
            }
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
                curr_ehh1_before_norm += del_update; // if(curr_ehh1_before_norm*1.0/normalizer_1 > p.EHH_CUTOFF){
            }else{
                curr_ehh0_before_norm += del_update;
            }
            curr_ehh_before_norm += del_update;
        }

        m.clear(); // CLEAR THE MAP //unordered_map<int, vector<int> >().swap(m);
        
        // if(curr_ehh1_before_norm*1.0/normalizer_1 > p.EHH_CUTOFF){
        // }

        // if(curr_ehh0_before_norm*1.0/normalizer_0 > p.EHH_CUTOFF){
        // }

        


        if(totgc == numHaps) {
            //std::cerr<<"Break reason for locus "<<locus<<":: ALL_UNIQUE."<<endl;
            break;
        }

        if(physicalDistance_from_core(i,locus, downstream) >= p.QWIN) {
            //std::cerr<<"Break reason for locus "<<locus<<":: MAX_EXTEND."<<endl;
            break;
        }
        // if(p.CALC_NSL && abs(i-locus) >= max_extend) {
        //     //std::cerr<<"Break reason for locus "<<locus<<":: MAX_EXTEND_NSL."<<endl;
        //     break; 
        // }


        
        double current_derived_ehh = curr_ehh1_before_norm*1.0/normalizer_1;
        double current_ancestral_ehh = curr_ehh0_before_norm*1.0/normalizer_0;
        double current_ehh = curr_ehh_before_norm*1.0/normalizer;
        output[i].ehh0 = current_ancestral_ehh;
        output[i].ehh1 = current_derived_ehh;
        output[i].ehh = current_ehh;
        output[i].pdist = downstream? -physicalDistance_from_core(i,locus, downstream): physicalDistance_from_core(i,locus, downstream);
        output[i].gdist = downstream? -geneticDistance_from_core(i,locus, downstream): geneticDistance_from_core(i,locus, downstream);
        output[i].print = true;

        // if(downstream){
        //     (*fout) << std::fixed <<   -(hm->mapData->mapEntries[locus].physicalPos -  hm->mapData->mapEntries[i].physicalPos)  << "\t"
        //     <<  -(hm->mapData->mapEntries[locus].geneticPos -  hm->mapData->mapEntries[i].geneticPos)<< "\t"
        //     << current_derived_ehh << "\t"
        //     << current_ancestral_ehh << "\t";
        // }else{
        //     (*fout) << std::fixed <<   hm->mapData->mapEntries[i].physicalPos -  hm->mapData->mapEntries[locus].physicalPos  << "\t"
        //     <<  hm->mapData->mapEntries[i].geneticPos -  hm->mapData->mapEntries[locus].geneticPos<< "\t"
        //     << current_derived_ehh << "\t"
        //     << current_ancestral_ehh << "\t";
        // }
        // (*fout) << current_ehh << endl;
    }

    delete[] group_count;
    delete[] group_id;
    delete[] isDerived;
}


/**
 * Calculate EHH in only one direction until cutoff is hit - upstream or downstream
*/
void EHH::calc_ehh_unidirection(int locus, unordered_map<int, vector<int> > & m, bool downstream){
    int numHaps = hm->hapData->nhaps;
    int numSnps = hm->mapData->nloci;
    bool unphased = p.UNPHASED;

    int total_iteration_of_m = 0;
    double ehh_before_norm = 0;
    double ehh0_before_norm = 0;
    double ehh1_before_norm = 0;

    bool gap_skip = false;

    uint64_t n_c0=0;
    uint64_t n_c1=0;

    uint64_t &ancestralCount = n_c0;
    uint64_t &derivedCount = n_c1;

    uint64_t core_n_c0=0;
    uint64_t core_n_c1=0;

    int group_count[numHaps];
    int group_id[numHaps];
    bool isDerived[numHaps];
    bool isAncestral[numHaps];

    //will be vectorized with compile time flags
    for(int i = 0; i<numHaps; i++){
        group_count[i] = 0;
        group_id[i] = 0;
        isDerived[i] = false;
        isAncestral[i] = false;
    }

    int totgc=0;
    vector<int> v = hm->hapData->hapEntries[locus].positions;
    
    if(v.size() == 0 or  v.size() == numHaps){
        std::cerr<<"ERROR: Monomorphic site should not exist";
        throw 0;
    }

    hm->hapData->hapEntries[locus].flipped = false;
    if(hm->hapData->hapEntries[locus].flipped){
        group_count[1] = v.size();
        group_count[0] = numHaps - v.size();
        n_c0 = v.size();
        n_c1 = numHaps - v.size();

        for (int set_bit_pos : v){
            isAncestral[set_bit_pos] = true;
            group_id[set_bit_pos] = 1;
        }
    }else{
        group_count[1] = v.size();
        group_count[0] = numHaps - v.size();
        n_c1 = v.size();
        n_c0 = numHaps - v.size();

        for (int set_bit_pos : v){
            isDerived[set_bit_pos] = true;
            group_id[set_bit_pos] = 1;
        }
    }        
    totgc+=2;
    ehh0_before_norm = twice_num_pair(n_c0);
    ehh1_before_norm = twice_num_pair(n_c1);
    ehh_before_norm = (twice_num_pair(n_c0)+twice_num_pair(n_c1));
    //twice_num_pair(n_c1+n_c0);

    
    int i = locus;  
    while(true){ // Upstream: for ( int i = locus+1; i<all_positions.size(); i++ )
        if(downstream){
            if (--i < 0) break;
            //if (hm.mentries[locus].phyPos - hm.mentries[i].phyPos > max_extend) break;
        }else{
            if (++i >= numSnps) break;
            //if (hm.mentries[i].phyPos -hm.mentries[locus].phyPos > max_extend) break;
        }
        
        
        //if(curr_ehh1_before_norm*1.0/n_c1_squared_minus < cutoff and curr_ehh0_before_norm*1.0/n_c0_squared_minus < cutoff){
        
        ///*
        if(ehh1_before_norm/twice_num_pair(n_c1) <= p.EHH_CUTOFF or ehh0_before_norm/twice_num_pair(n_c0)  <= p.EHH_CUTOFF){   // or cutoff, change for benchmarking against hapbin
            //std::cout<<"breaking"<<endl;
            break;
        }
        //*/
        double distance;
        
        if(downstream){
            distance = hm->mapData->mapEntries[i+1].physicalPos - hm->mapData->mapEntries[i].physicalPos;
        }else{
            distance = hm->mapData->mapEntries[i].physicalPos - hm->mapData->mapEntries[i-1].physicalPos;
        }

        // this should not happen as we already did integrity check previously
        if (distance < 0)
        {
            std::cerr << "ERROR: physical position not in ascending order.\n"; 
            throw 0;
        }
        
        
        if(downstream){
            v = hm->hapData->hapEntries[i+1].xors;
        }else{
            v = hm->hapData->hapEntries[i].xors;
        }

        if(hm->hapData->hapEntries[i].positions.size() < v.size() && i!=numHaps-1 ){ //  dont do in boundary
            v = hm->hapData->hapEntries[i].positions;
        }
        
        for (const int& set_bit_pos : v){
            int old_group_id = group_id[set_bit_pos];
            m[old_group_id].push_back(set_bit_pos);
        }

        for (const auto &ele : m) {
            
            int old_group_id = ele.first;
            int newgroup_size = ele.second.size() ;
                            
            total_iteration_of_m += newgroup_size;
                            
            if(group_count[old_group_id] == newgroup_size || newgroup_size == 0){
                continue;
            }

            for(int v: ele.second){
                group_id[v] = totgc;
            }
            
            double del_update = -twice_num_pair(group_count[old_group_id]) + twice_num_pair(newgroup_size) + twice_num_pair(group_count[old_group_id] - newgroup_size);
            if(p.ALT){
                del_update = -square_alt(group_count[old_group_id]) +   square_alt(newgroup_size) + square_alt(group_count[old_group_id] - newgroup_size);
            }
            
            group_count[old_group_id] -= newgroup_size;
            group_count[totgc] += newgroup_size;
            
            totgc+=1;
            
            ehh_before_norm += del_update;

            //bool isDerivedGroup =  (!hm->hapData->hapEntries[locus].flipped && isDerived[ele.second[0]]) || (hm->hapData->hapEntries[locus].flipped && !isAncestral[ele.second[0]]); // just check first element to know if it is derived. 
                bool isDerivedGroup = isDerived[ele.second[0]];
            if(isDerivedGroup) // if the core locus for this chr has 1 (derived), then update ehh1, otherwise ehh0
            {
                ehh1_before_norm += del_update;
            }else{
                ehh0_before_norm += del_update;
            }
        }
        m.clear();
        

        //printing 

        double current_derived_ehh; 
        double current_ancestral_ehh;
        double current_ehh;
        
        if(p.ALT){
            current_ehh = ehh_before_norm*1.0/square_alt(n_c1+n_c0);
            current_derived_ehh = ehh1_before_norm*1.0/square_alt(n_c1);
            current_ancestral_ehh  = ehh0_before_norm*1.0/square_alt(n_c0);
            
        }else{
            current_ehh = ehh_before_norm*1.0/twice_num_pair(n_c1+n_c0);
            current_derived_ehh = ehh1_before_norm*1.0/twice_num_pair(n_c1);
            current_ancestral_ehh  = ehh0_before_norm*1.0/twice_num_pair(n_c0);
            
        }

        if(downstream){
            (*fout) << std::fixed <<   -int(hm->mapData->mapEntries[locus].physicalPos -  hm->mapData->mapEntries[i].physicalPos)  << "\t"
            <<  -(hm->mapData->mapEntries[locus].geneticPos -  hm->mapData->mapEntries[i].geneticPos)<< "\t"
            << current_derived_ehh << "\t"
            << current_ancestral_ehh << "\t";
        }else{
            (*fout) << std::fixed <<   hm->mapData->mapEntries[i].physicalPos -  hm->mapData->mapEntries[locus].physicalPos  << "\t"
            <<  hm->mapData->mapEntries[i].geneticPos -  hm->mapData->mapEntries[locus].geneticPos<< "\t"
            << current_derived_ehh << "\t"
            << current_ancestral_ehh << "\t";
        }

        /*
        cout<<"Iter "<<i-locus<<": EHH1["<<locus<<","<<i<<"]=";

        for (int x = 0 ; x < totgc;  x++){
            cout<<group_count[x]<<"("<<x<<")";
        }

        //print all elements of vector v
        for (int x = 0 ; x < v.size();  x++){
            cout<<v[x]<<" ";
        }
        
        cout<<endl;
        */

            
        if(unphased){
            //(*fout) << current_notAncestral_ehh << "\t"
            //        << current_notDerived_ehh << "\t";
        }
        (*fout) << current_ehh << endl;
    }
}


/**
 * @brief Calculate the EHH for a single locus
 * @param query The query locus name
*/
void EHH::calc_single_ehh(string query){
    int numSnps = hm->mapData->nloci;
    int numHaps = hm->hapData->nhaps;

    int locus = hm->queryFound(query);

    if(p.UNPHASED){
        throw "Unphased EHH not implemented yet";
    }else{
        calc_ehh_unidirection(locus, true); // downstream
        calc_ehh_unidirection(locus, false); // upstream
        //now print

        (*fout) << std::fixed <<   "pdist\tgdist\tderEHH\tancEHH\tEHH" << "\n";
        for (int i = 0; i < numSnps; i++){
            if(!output[i].print){
                continue;
            }
            (*fout) << std::fixed <<   output[i].pdist  << "\t"
            <<  output[i].gdist << "\t"
            << output[i].ehh1 << "\t"
            << output[i].ehh0  << "\t"
            << output[i].ehh  << "";
            (*fout) << endl;

            // if(downstream){
            //     (*fout) << std::fixed <<   -int(hm->mapData->mapEntries[locus].physicalPos -  hm->mapData->mapEntries[i].physicalPos)  << "\t"
            //     <<  -(hm->mapData->mapEntries[locus].geneticPos -  hm->mapData->mapEntries[i].geneticPos)<< "\t"
            //     << current_derived_ehh << "\t"
            //     << current_ancestral_ehh << "\t";
            // }else{
            //     (*fout) << std::fixed <<   hm->mapData->mapEntries[i].physicalPos -  hm->mapData->mapEntries[locus].physicalPos  << "\t"
            //     <<  hm->mapData->mapEntries[i].geneticPos -  hm->mapData->mapEntries[locus].geneticPos<< "\t"
            //     << current_derived_ehh << "\t"
            //     << current_ancestral_ehh << "\t";
            // }
        }
    }
}

// /**
//  * @brief Calculate the EHH for a single locus
//  * @param query The query locus phys pos
// */
// void EHH::calc_single_ehh(int query){
//     int numSnps = hm->mapData->nloci;
//     int numHaps = hm->hapData->nhaps;

//     //TODO
//     int locus = hm.queryFound(query);
    
//     unordered_map<int, vector<int> > m;
//     calc_ehh_unidirection(locus, m, false); // upstream
//     calc_ehh_unidirection(locus, m, true); // downstream
    
// }