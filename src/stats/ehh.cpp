#include "ehh.h"
#include <cassert>

#define UNDEFINED_EHH -1

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

        if(downstream)
        (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName << " is monomorphic. Undefined EHH is reported as "<<UNDEFINED_EHH <<"."<<endl;
        // (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
        //         << " (number " << locus << ") is monomorphic. Skipping calculation at this locus. "<<endl;
        // output[locus].print = false;
        // return;
    }else if (n_c1==numHaps){ // all 1s

        if(downstream)
        (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName << " is monomorphic. Undefined EHH is reported as "<<UNDEFINED_EHH <<"."<<endl;

        group_count[0] = numHaps;
        totgc+=1;
        
        ACTION_ON_ALL_SET_BITS(get_all_1s(locus), {
            isDerived[set_bit_pos] = true;
        });
        
        curr_ehh1_before_norm = normalizer_1;
        // (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
        //         << " (number " << locus << ") is monomorphic. Skipping calculation at this locus. "<<endl;
        // output[locus].print = false;
        // return;
    }else{  //so both n_c1 and n_c0 is non-0
        group_count[1] = n_c1;
        group_count[0] = n_c0;
        
        ACTION_ON_ALL_SET_BITS(get_all_1s(locus), {
            isDerived[set_bit_pos] = true;
            group_id[set_bit_pos] = 1;
        });
        
        totgc+=2;
        curr_ehh0_before_norm = normalizer_0;
        curr_ehh1_before_norm = normalizer_1;
    }
    curr_ehh_before_norm = normalizer_0 + normalizer_1;


    output[locus].ehh0 = (normalizer_0==0)? UNDEFINED_EHH : curr_ehh0_before_norm*1.0/normalizer_0;
    output[locus].ehh1 = (normalizer_1==0)? UNDEFINED_EHH : curr_ehh1_before_norm*1.0/normalizer_1;
    output[locus].ehh = curr_ehh_before_norm*1.0/normalizer;
    output[locus].gdist = 0;
    output[locus].pdist = 0;
    output[locus].print = true;
    
    //cout<<locus<< " "<<output[locus].pdist<<" "<<output[locus].ehh0<<" "<<output[locus].ehh1<<" "<<output[locus].ehh<<endl;

    int i = locus;  
    while(true){ // Upstream: for ( int i = locus+1; i<all_positions.size(); i++ )
        //output[i].print = false;

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
                //std::cout<<"Break reason for locus "<<locus<<":: REACHED_LEFT_EDGE."<<endl;
                break;
            }
        }else{
            if (++i >= numSnps) {
                //std::cout<<"Break reason for locus "<<locus<<":: REACHED_RIGHT_EDGE."<<endl;
                break;
            }
        }

        double distance =  geneticDistance(i, downstream);
        if (distance < 0) // this should not be needed as we already did integrity check previously
        {
            std::cerr << "ERROR: physical position not in ascending order.\n"; 
            exit(EXIT_FAILURE);
        }
        // if(p.CALC_NSL){
        //     distance = 1;
        // }
        double scale = double(p.SCALE_PARAMETER) / double(physicalDistance(i, downstream) );
        if (scale > 1) scale = 1;
        distance *= scale;

        if(hm->hapData->get_n_c0(i) == 0 or hm->hapData->get_n_c1(i) == 0 ){ // monomorphic check, do not compute unnecessarily
        //if(false){
        // pthread_mutex_lock(&mutex_log);
            // std::cerr<<"ERROR: Monomorphic site should not exist."<<endl;
            // std::cerr<< hm->hapData->get_n_c0(i) <<" n_c0 at locus "<< i <<endl; 
            // std::cerr<< hm->hapData->get_n_c1(i) <<" n_c1 at locus "<< i<< endl; 
            // pthread_mutex_unlock(&mutex_log);
            // throw 0;
            
            //if you wish to continue anyway
            double current_derived_ehh = normalizer_1==0? UNDEFINED_EHH: curr_ehh1_before_norm*1.0/normalizer_1;
            double current_ancestral_ehh = normalizer_0==0? UNDEFINED_EHH: curr_ehh0_before_norm*1.0/normalizer_0;
            double current_ehh = curr_ehh_before_norm*1.0/normalizer;

            output[i].ehh0 = current_ancestral_ehh;
            output[i].ehh1 = current_derived_ehh;
            output[i].ehh = current_ehh;
            output[i].pdist = downstream? -physicalDistance_from_core(i,locus, downstream): physicalDistance_from_core(i,locus, downstream);
            output[i].gdist = downstream? -geneticDistance_from_core(i,locus, downstream): geneticDistance_from_core(i,locus, downstream);
            output[i].print = true;
            //cout<<i<< " "<<output[i].pdist<<" "<<output[i].ehh0<<" "<<output[i].ehh1<<" "<<output[i].ehh<<endl;
        }else{
            // PHASE 2: GET MAP
            std::unique_ptr<std::unordered_map<int, std::vector<int>>> mp(new std::unordered_map<int, std::vector<int>>());
            unordered_map<int, vector<int> >& m = (* mp);
            ACTION_ON_ALL_SET_BITS(get_all_1s(i, downstream), {
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
            double current_derived_ehh = normalizer_1==0?  UNDEFINED_EHH: curr_ehh1_before_norm*1.0/normalizer_1;
            double current_ancestral_ehh = normalizer_0==0? UNDEFINED_EHH: curr_ehh0_before_norm*1.0/normalizer_0;
            double current_ehh = curr_ehh_before_norm*1.0/normalizer;
            output[i].ehh0 = current_ancestral_ehh;
            output[i].ehh1 = current_derived_ehh;
            output[i].ehh = current_ehh;
            output[i].pdist = downstream? -physicalDistance_from_core(i,locus, downstream): physicalDistance_from_core(i,locus, downstream);
            output[i].gdist = downstream? -geneticDistance_from_core(i,locus, downstream): geneticDistance_from_core(i,locus, downstream);
            output[i].print = true;
            //cout<<i<< " "<<output[i].pdist<<" "<<output[i].ehh0<<" "<<output[i].ehh1<<" "<<output[i].ehh<<endl;

        }

        // if(totgc == numHaps) {
        //     //std::cerr<<"Break reason for locus "<<locus<<":: ALL_UNIQUE."<<endl;
        //     break;
        // }

        if(physicalDistance_from_core(i,locus, downstream) > p.QWIN) {
            //std::cout<<"Break reason for locus "<<locus<<":: MAX_EXTEND." << physicalDistance_from_core(i,locus, downstream)  <<endl;
            output[i].print = false;
            break;
        }

    }

    delete[] group_count;
    delete[] group_id;
    delete[] isDerived;
}


/**
 * @brief Calculate the EHH for a single locus
 * @param query The query locus name
*/
void EHH::main(string query){
    init_global_fout("ehh");
    init_output_and_querymap();

    int numSnps = hm->mapData->nloci;
    int numHaps = hm->hapData->nhaps;

    int locus = this->queryFound(query);
    if(locus == -1){
        return;
    }

    if(p.UNPHASED){
        calc_ehh_unidirection_unphased(locus, true); // downstream
        calc_ehh_unidirection_unphased(locus, false); // upstream

        (*fout) << std::fixed <<   "pdist\tgdist\tehh2\tehh0\tcehh\tehh" << "\n";
        for (int i = 0; i < numSnps; i++){
            if(!output[i].print){
                continue;
            }
            (*fout) << std::fixed <<   output[i].pdist  << "\t"
            <<  output[i].gdist << "\t"
            << output[i].ehh2 << "\t"
            << output[i].ehh0  << "\t"
            << output[i].cehh  << "\t"
            << output[i].ehh  << "";
            (*fout) << endl;
        }
        
    }else{
        for (int i = 0; i < numSnps; i++){
            (output[i].print = false);
        }

        //vector<EHHOutput> output(numSnps);
        calc_ehh_unidirection(locus, true); // downstream
        // for (int i = 0; i < numSnps; i++){
        //     if(output[i].print == false){
        //         continue;
        //     }
        //     // (cout) << std::fixed <<   output[i].pdist  << "\t"
        //     // <<  output[i].gdist << "\t"
        //     // << output[i].ehh1 << "\t"
        //     // << output[i].ehh0  << "\t"
        //     // << output[i].ehh  << "";
        //     // (cout) << endl;
        //     (*fout) << std::fixed <<   output[i].pdist  << "\t"
        //     <<  output[i].gdist << "\t"
        //     << output[i].ehh1 << "\t"
        //     << output[i].ehh0  << "\t"
        //     << output[i].ehh  << "";
        //     (*fout) << endl;
        // }

        calc_ehh_unidirection(locus, false); // upstream

        (*fout) << std::fixed <<   "pdist\tgdist\tderEHH\tancEHH\tEHH" << "\n";
        for (int i = 0; i < numSnps; i++){
            if(output[i].print == false){
                continue;
            }
            // (cout) << std::fixed <<   output[i].pdist  << "\t"
            // <<  output[i].gdist << "\t"
            // << output[i].ehh1 << "\t"
            // << output[i].ehh0  << "\t"
            // << output[i].ehh  << "";
            // (cout) << endl;
            (*fout) << std::fixed <<   output[i].pdist  << "\t"
            <<  output[i].gdist << "\t"
            << output[i].ehh1 << "\t"
            << output[i].ehh0  << "\t"
            << output[i].ehh  << "";
            (*fout) << endl;
        }
    }
}


void EHH::calc_ehh_unidirection_unphased(int locus, bool downstream){

    std::unique_ptr<std::unordered_map<int, std::vector<int>>> mp(new std::unordered_map<int, std::vector<int>>());
    unordered_map<int, vector<int> >& m = (* mp);

    int numSnps = hm->hapData->nloci;
    int numHaps = hm->hapData->nhaps;

    double curr_ehh_before_norm[3];
    double curr_cehh_before_norm[3];
    double ehh_before_norm = 0;

    int n_c[3] = {0,0,0};

    int* group_count = new int[numHaps];
    int* group_id =  new int[numHaps];
    int* group_core =  new int[numHaps];
    bool* is1 =  new bool[numHaps];
    bool* is2 =  new bool[numHaps];
    
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
        cerr<<"ERROR: n_c1 + n_c2 + n_c0 != numHaps"<<endl;
        exit(2);
    }

    string orderStr = getOrder(n_c[2], n_c[1], n_c[0]);

    double normalizer[3];
    double normalizer_not[3];
    double normalizer_just = twice_num_pair_or_square(numHaps, p.ALT);
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

    ACTION_ON_ALL_SET_BITS(get_all_1s(locus), {
        is1[set_bit_pos] = true;
        group_id[set_bit_pos] = pos_of_012[1];   
    });

    ACTION_ON_ALL_SET_BITS(get_all_2s(locus), {
        is2[set_bit_pos] = true;
        group_id[set_bit_pos] = pos_of_012[2]; 
    });
   
    if(group_count[0] == numHaps){ //monomorphic site
        totgc+=1;
    }else if(group_count[2] == 0 ){ //second clause is redundant
        if(group_count[0] + group_count[1] != numHaps){
            cerr<<"ERROR: gc2==0 locus"<<locus<<endl;
            exit(2);
        }
        totgc+=2;
    }else{
        totgc+=3;
    }

    if(n_c[1] == numHaps || n_c[0] == numHaps || n_c[2] == numHaps){ 
        (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
                << " (number " << locus << ") is monomorphic. Skipping calculation at this locus. "
                << "het: " <<  n_c[1] << " hom0: " << n_c[0] << " hom1: " << n_c[2] << ".\n";
        return;
    }

    double freqHetGT = n_c[1]*1.0/numHaps;
    if (  freqHetGT > 1-p.MAF ) 
    {
        (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
                << " (number " << locus << ") has too many hets. Skipping calculation at this locus. "
                << "het: " <<  n_c[1] << " hom0: " << n_c[0] << " hom1: " << n_c[2] << ".\n";
        return;
    }

    output[locus].ehh0 = curr_ehh_before_norm[0]*1.0/normalizer[0];
    output[locus].ehh2 = curr_ehh_before_norm[2]*1.0/normalizer[2];
    output[locus].ehh = ehh_before_norm*1.0/normalizer_just;
    output[locus].gdist = 0;
    output[locus].pdist = 0;
    output[locus].print = true;
    

    int i = locus;  // locus == core_locus
    while(true){ // Upstream: for ( int i = locus+1; i<all_positions.size(); i++ )
    
        // if(p.CALC_IHS && !p.CALC_NSL){
        //     if(curr_ehh_before_norm[2]*1.0/normalizer[2] <= p.EHH_CUTOFF and curr_ehh_before_norm[0]*1.0/normalizer[0]  <= p.EHH_CUTOFF){   // or cutoff, change for benchmarking against hapbin
        //         //std::cerr<<"Break reason for locus "<<locus<<":: EHH_CUTOFF."<<endl;
        //         break;
        //     }
        // }

        bool edgeBreak = false;
        edgeBreak = nextLocOutOfBounds(i, downstream);
        if(edgeBreak) {
            // (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << p.EHH_CUTOFF
            //         << ". position: "<< hm->mapData->mapEntries[locus].physicalPos << " id: " << hm->mapData->mapEntries[locus].locusName << endl;
            break;
        }
        i = (downstream) ? i-1 : i+1;
        
        double distance =  geneticDistance(i, downstream);
        if (distance < 0) // this should not be needed as we already did integrity check previously
        {
            std::cerr << "ERROR: physical position not in ascending order.\n"; 
            exit(2);
        }

        // if (physicalDistance(i,downstream) > p.MAX_GAP)
        // {
        //     (*flog) << "WARNING: Reached a gap of " << physicalDistance(i,downstream)
        //             << "bp > " << p.MAX_GAP << "bp. Skipping calculation at position " << hm->mapData->mapEntries[locus].physicalPos << " id: " << hm->mapData->mapEntries[locus].locusName << "\n";
        //     return;
        // }


       if(hm->hapData->get_n_c0(i) == numHaps or hm->hapData->get_n_c1(i) == numHaps or hm->hapData->get_n_c2(i) == numHaps){ // monomorphic check, do not compute unnecessarily
            // pthread_mutex_lock(&mutex_log);
            // std::cerr<<"ERROR: Monomorphic site should not exist."<<endl;
            // std::cerr<< hm->hapData->get_n_c0(i) <<" n_c0 at locus "<< i <<endl; 
            // std::cerr<< hm->hapData->get_n_c1(i) <<" n_c1 at locus "<< i<< endl; 
            // pthread_mutex_unlock(&mutex_log);
            // throw 0;
            
            //if you wish to continue anyway

            continue; 
        }

        
        ACTION_ON_ALL_SET_BITS(get_all_2s(i), {
            int old_group_id = group_id[set_bit_pos];
            m[old_group_id].push_back(set_bit_pos);
        });
        updateEHH_from_split_unphased(m, group_count, group_id, totgc, curr_ehh_before_norm, curr_cehh_before_norm, is1, is2, group_core);
        m.clear();

        ACTION_ON_ALL_SET_BITS(get_all_1s(i), {
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
        double ehh_before_norm = 0;
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
            ehh_before_norm += gcsquare;
        }

        if(double(curr_ehh_before_norm[0]*1.0/normalizer[0]) > p.EHH_CUTOFF){
            curr_ehh_before_norm[0] = ehh0_before_norm;
        }

        if(double(curr_ehh_before_norm[2]*1.0/normalizer[2]) > p.EHH_CUTOFF ){
            curr_ehh_before_norm[2] = ehh2_before_norm;
        }
            
        if(true){  //this is how it's in selscan, does not depend on cutoff
            curr_cehh_before_norm[0] = cehh0_before_norm;
            curr_cehh_before_norm[2] = cehh2_before_norm;
        }

        //if(totgc == numHaps) {
            // std::cerr<<"DEBUG::: Break reason for locus "<<locus<<":: ALL_UNIQUE."<<endl;
            // break;
        //}

        if(physicalDistance_from_core(i,locus, downstream) >= p.QWIN) {
            //std::cerr<<"Break reason for locus "<<locus<<":: MAX_EXTEND."<<endl;
           break;
        }

        output[i].ehh0 = curr_ehh_before_norm[0]*1.0/normalizer[0];
        output[i].ehh2 = curr_ehh_before_norm[2]*1.0/normalizer[2];
        output[i].ehh = ehh_before_norm*1.0/normalizer_just;
        output[i].pdist = downstream? -physicalDistance_from_core(i,locus, downstream): physicalDistance_from_core(i,locus, downstream);
        output[i].gdist = downstream? -geneticDistance_from_core(i,locus, downstream): geneticDistance_from_core(i,locus, downstream);
        output[i].print = true;

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
}