#include "ihs.h"
#include <set>
#include <algorithm>
#include <cassert>

#include <pthread.h>


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



pthread_mutex_t IHS::mutex_log = PTHREAD_MUTEX_INITIALIZER;

void IHS::updateEHH_from_split_unphased( map<int, vector<int> >& m, int* group_count, int* group_id, int& totgc, uint64_t* ehh_before_norm, uint64_t* cehh_before_norm, bool* is1, bool* is2){
    for (const auto &ele : m) {
        int old_group_id = ele.first;
        int newgroup_size = ele.second.size() ;
        if(group_count[old_group_id] == newgroup_size || newgroup_size == 0){
            continue;
        }
        for(int v: ele.second){
            group_id[v] = totgc;
        }
        
        double del_update = -twice_num_pair_or_square(group_count[old_group_id], p.ALT) + twice_num_pair_or_square(newgroup_size, p.ALT) + twice_num_pair_or_square(group_count[old_group_id] - newgroup_size, p.ALT);
        // if(p.ALT){
        //     del_update = -square_alt(group_count[old_group_id]) +   square_alt(newgroup_size) + square_alt(group_count[old_group_id] - newgroup_size);
        // }
        
        group_count[old_group_id] -= newgroup_size;
        //ifgroup_count[old_group_id] == 0 , dont inc totgc
        group_count[totgc] += newgroup_size;
        
        totgc+=1;
        
        //TODO
        if(is1[ele.second[0]]) // if the core locus for this chr has 1 (derived), then update ehh1, otherwise ehh0
        {
            ehh_before_norm[1] += del_update;

            cehh_before_norm[2] += del_update;
            cehh_before_norm[0] += del_update;
        }else if(is2[ele.second[0]]) {
            ehh_before_norm[2] += del_update;

            cehh_before_norm[1] += del_update;
            cehh_before_norm[0] += del_update;
        }else{
            ehh_before_norm[0] += del_update;

            cehh_before_norm[2] += del_update;
            cehh_before_norm[1] += del_update;
        }
    }
}

void IHS::calc_ehh_unidirection_unphased(int locus, bool downstream){
    int numSnps = hm->hapData->nloci;
    int numHaps = hm->hapData->nhaps;

    uint64_t prev_ehh_before_norm[3] = {0,0,0};
    uint64_t curr_ehh_before_norm[3] = {0,0,0};

    uint64_t prev_cehh_before_norm[3] = {0,0,0};
    uint64_t curr_cehh_before_norm[3] = {0,0,0};

    int n_c[3] = {0,0,0};

    bool gap_skip = false;

    int group_count[numHaps];
    int group_id[numHaps];
    bool is1[numHaps];
    bool is2[numHaps];

    double* iHH[3] = {iHH0, iHH1, iHH2};
    double* ciHH[3] = {ciHH0, ciHH1, ciHH2};

    //will be vectorized with compile time flags
    for(int i = 0; i<numHaps; i++){
        group_count[i] = 0;
        is1[i] = false;
        is2[i] = false;
    }

    int totgc=0;

    vector<unsigned int> v1 = hm->hapData->hapEntries[locus].positions;
    vector<unsigned int> v2 = hm->hapData->hapEntries[locus].positions2;

    vector<unsigned int> g0 = hm->hapData->hapEntries[locus].g[0];
    vector<unsigned int> g1 = hm->hapData->hapEntries[locus].g[1];
    vector<unsigned int> g2 = hm->hapData->hapEntries[locus].g[2];
  

    n_c[0] = numHaps - v1.size() - v2.size();
    n_c[2] = v2.size();
    n_c[1] = v1.size();

    
    /*
    n_c[0] = numHaps - v1.size();
    n_c[2] = 0;
    n_c[1] = 0;
    for(const int& set_bit_pos : v1){
        if(set_bit_pos%2){
            n_c[2]++;
        }else{
            n_c[1]++;
        }
    }
    */



    if(n_c[1] + n_c[2] + n_c[0] != numHaps){
        cerr<<"ERROR: n_c1 + n_c2 + n_c0 != numHaps"<<endl;
        exit(1);
    }

    // n_c[2] = v2.size();
    // n_c[1] = v1.size();
    // n_c[0] = numHaps - n_c[1] - n_c[2];
    
    
    string orderStr = getOrder(n_c[2], n_c[1], n_c[0]);
    for(int i = 0; i<3; i++){
        curr_ehh_before_norm[i] = twice_num_pair_or_square(n_c[i], p.ALT);
        curr_cehh_before_norm[i] = twice_num_pair_or_square(numHaps -  n_c[i], p.ALT);
        group_count[i] = n_c[orderStr[i]-'0'];
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

    /*
    for(const int& set_bit_pos : v1){
        if(set_bit_pos%2){
            is2[set_bit_pos/2] = true;
            group_id[set_bit_pos/2] = pos_of_012[2];
        }else{
            is1[set_bit_pos/2] = true;
            group_id[set_bit_pos/2] = pos_of_012[1];
        }
    }
    */

   for(const int& set_bit_pos : v2){
        is2[set_bit_pos] = true;
        group_id[set_bit_pos] = pos_of_012[2];
    }
    for(const int& set_bit_pos : v1){
        is1[set_bit_pos] = true;
        group_id[set_bit_pos] = pos_of_012[1];
    }

    //TODO
    // if(group_count[0] == numHaps){ //monomorphic site
    //     totgc+=1;
    // }else if(group_count[2] == 0 ){ //second clause is redundant
    //     if(group_count[1] + group_count[2] != numHaps){
    //         cerr<<"ERROR: gc2==0"<<endl;
    //         exit(1);
    //     }
    //     totgc+=2;
    // }else{
    //     totgc+=3;
    // }

    if(group_count[0] == numHaps){ //monomorphic site
        totgc+=1;
    }else if(group_count[2] == 0 ){ //second clause is redundant
        if(group_count[0] + group_count[1] != numHaps){
            cerr<<"ERROR: gc2==0 locus"<<locus<<endl;
            exit(1);
        }
        totgc+=2;
    }else{
        totgc+=3;
    }

    if(downstream){
        for(int i=0; i<3; i++){
            if(twice_num_pair_or_square(n_c[i], p.ALT)!=0){
                iHH[i][locus] += (curr_ehh_before_norm[i] + prev_ehh_before_norm[i]) * 0.5 / twice_num_pair_or_square(n_c[i], p.ALT);
            }

            if(twice_num_pair_or_square(numHaps - n_c[i], p.ALT)!=0){
                ciHH[i][locus] += (curr_cehh_before_norm[i] + prev_cehh_before_norm[i])  * 0.5 / twice_num_pair_or_square(numHaps - n_c[i], p.ALT);
                
            }
        }
    }
    
    for(int i=0; i<3; i++){
        prev_ehh_before_norm[i] = curr_ehh_before_norm[i];
        prev_cehh_before_norm[i] =  curr_cehh_before_norm[i];
    }
    

    int i = locus;  
    while(true){ // Upstream: for ( int i = locus+1; i<all_positions.size(); i++ )
        if(downstream){
            if (--i < 0) break;
            //if (hm.mentries[locus].phyPos - hm.mentries[i].phyPos > max_extend) break;
        }else{
            if (++i >= numSnps) break;
            //if (hm.mentries[i].phyPos -hm.mentries[locus].phyPos > max_extend) break;
        }
        
        
        if(curr_ehh_before_norm[2]*1.0/twice_num_pair_or_square(n_c[2], p.ALT) <= p.EHH_CUTOFF and curr_ehh_before_norm[0]*1.0/twice_num_pair_or_square(n_c[0], p.ALT)  <= p.EHH_CUTOFF){   // or cutoff, change for benchmarking against hapbin
            //std::cerr<<"Break reason for locus "<<locus<<":: EHH_CUTOFF."<<endl;
            break;
        }

        
        //TODO
        // if (current_derived_ehh <=  EHH_CUTOFF) dont add to ihh calculation, only do EHH ancestral 
        //     {
        //     }


        // double distance;
        
        // if(downstream){
        //     distance = hm->mapData->mapEntries[i+1].physicalPos - hm->mapData->mapEntries[i].physicalPos;
            
        // }else{
        //     distance = hm->mapData->mapEntries[i].physicalPos - hm->mapData->mapEntries[i-1].physicalPos;
        // }
        // double scale = double(p.SCALE_PARAMETER) / distance;
        // if (scale > 1) scale = 1;
        // //distance = 1; // for testing
        
        // if(p.CALC_NSL){
        //     distance = 1;
        //     scale = 1;
        // }



        double scale, distance;
        if(p.CALC_NSL){
            scale = double(p.SCALE_PARAMETER) / geneticDistance(i, downstream);  
            distance = geneticDistance(i, downstream);
            distance = 1;

        }else{
            scale = double(p.SCALE_PARAMETER) / physicalDistance(i, downstream);   
            distance = physicalDistance(i, downstream);  
        }
        if(distance > p.SCALE_PARAMETER){
            distance /= p.SCALE_PARAMETER;
        }
        
        if(scale > 1) scale = 1;

        
        // this should not happen as we already did integrity check previously
        if (distance < 0)
        {
            std::cerr << "ERROR: physical position not in ascending order.\n"; 
            throw 0;
        }
        //if (distance >= p.MAX_EXTEND) break;
        
        vector<unsigned int> xors;
        if(downstream){
            g0 = hm->hapData->hapEntries[i+1].g[0];
            g1 = hm->hapData->hapEntries[i+1].g[1];
            g2 = hm->hapData->hapEntries[i+1].g[2];
        }else{
            g0 = hm->hapData->hapEntries[i].g[0];
            g1 = hm->hapData->hapEntries[i].g[1];
            g2 = hm->hapData->hapEntries[i].g[2];
        }

        // if(distance> max_gap){
        //     gap_skip = true;
        //     break;
        // }


        // if(distance > p.SCALE_PARAMETER){
        //     distance /= p.SCALE_PARAMETER;
        // }

        
        if(n_c[0] == numHaps or n_c[1] == numHaps or n_c[2] == numHaps){
            cerr<<"WARNING: Monomorphic site."<<endl;
        
            for(int i=0; i<3; i++){
                if(twice_num_pair_or_square(n_c[i], p.ALT)!=0){
                    iHH[i][locus] += (curr_ehh_before_norm[i] + prev_ehh_before_norm[i])  * distance * 0.5 / twice_num_pair_or_square(n_c[i], p.ALT);
                }
                if(twice_num_pair_or_square(numHaps - n_c[i], p.ALT)!=0){
                    //if(cehh_before_norm[i]*1.0/twice_num_pair(numHaps - n_c[i]) > p.EHH_CUTOFF){
                        ciHH[i][locus] += (curr_cehh_before_norm[i] + prev_cehh_before_norm[i])  * distance * 0.5 / twice_num_pair_or_square(numHaps - n_c[i], p.ALT);
                    //}
                }
            }

            
            
            continue;
        }
        //cout << locus<<" " << iHH2[locus]<<" "<<iHH1[locus]<<" "<<iHH0[locus]<<endl;

        map<int, vector<int> > m;
        for (const unsigned int& set_bit_pos : g0){
            int old_group_id = group_id[set_bit_pos];
            m[old_group_id].push_back(set_bit_pos);         
        }
        updateEHH_from_split_unphased(m, group_count, group_id, totgc, curr_ehh_before_norm, curr_cehh_before_norm, is1, is2);
        m.clear();

        for (const unsigned int& set_bit_pos : g1){
            int old_group_id = group_id[set_bit_pos];
            m[old_group_id].push_back(set_bit_pos);         
        }
        updateEHH_from_split_unphased(m, group_count, group_id, totgc, curr_ehh_before_norm, curr_cehh_before_norm, is1, is2);
        m.clear();


        for (const unsigned int& set_bit_pos : g2){
            int old_group_id = group_id[set_bit_pos];
            m[old_group_id].push_back(set_bit_pos);         
        }
        updateEHH_from_split_unphased(m, group_count, group_id, totgc, curr_ehh_before_norm, curr_cehh_before_norm, is1, is2);
        m.clear();
        

        // if(twice_num_pair(n_c1)!=0){

        //     if(ehh1_before_norm*1.0/twice_num_pair(n_c1) > p.EHH_CUTOFF){
               
        //         iHH1[locus] += (curr_ehh1_before_norm + ehh1_before_norm) * distance * 0.5 / twice_num_pair(n_c1);
        //     }
            
        //     //cout<<"Summing "<<1.0*curr_ehh1_before_norm/n_c1_squared_minus<<" "<<1.0*ehh1_before_norm/n_c1_squared_minus<<endl;
        // }

        // if(twice_num_pair(n_c0)!=0){
        //     if(ehh0_before_norm*1.0/twice_num_pair(n_c0) > p.EHH_CUTOFF){
        //         iHH0[locus] += (curr_ehh0_before_norm + ehh0_before_norm) * distance * 0.5 / twice_num_pair(n_c0);

        //     }
        // }

        for(int i=0; i<3; i++){
            if(twice_num_pair_or_square(n_c[i], p.ALT)!=0){
                if(curr_ehh_before_norm[i]*1.0/twice_num_pair_or_square(n_c[i], p.ALT) > p.EHH_CUTOFF){
                    // this was working iHH[i][locus] += (curr_ehh_before_norm[i] + prev_ehh_before_norm[i]) * scale * distance* 0.5 / twice_num_pair(n_c[i]);
                    iHH[i][locus] += (curr_ehh_before_norm[i] + prev_ehh_before_norm[i])  * distance * 0.5 / twice_num_pair_or_square(n_c[i], p.ALT);

                }
            }


            if(twice_num_pair_or_square(numHaps - n_c[i], p.ALT)!=0){
                //this is how it's in selscan //if(cehh_before_norm[i]*1.0/twice_num_pair(numHaps - n_c[i]) > p.EHH_CUTOFF){
                    ciHH[i][locus] += (curr_cehh_before_norm[i] + prev_cehh_before_norm[i])  * distance * 0.5 / twice_num_pair_or_square(n_c[i], p.ALT);
                //}
            }

            prev_ehh_before_norm[i] = curr_ehh_before_norm[i];
            prev_cehh_before_norm[i] = curr_cehh_before_norm[i];
        }
        //cout << locus<<" " << iHH2[locus]<<" "<<iHH1[locus]<<" "<<iHH0[locus]<<" " << ciHH2[locus]<<" " << ciHH1[locus]<<" " << ciHH0[locus]<<endl;


        // pthread_mutex_lock(&mutex_log);
        // for(int i =0; i<numHaps; i++){
        //     if(group_count[i] == 0){
        //          cout<<i<<" ";
        //     }
           
        // }
        // pthread_mutex_unlock(&mutex_log);
        // cout<<endl;
        // cout<<endl;
        // cout<<endl;


        // if (nextLocus < 0 || nextLocus >= hm->mapData->nloci)
        // {
        //     pthread_mutex_lock(&mutex_log);
        //     (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << p.EHH_CUTOFF
        //             << ". ";
        //     if (!p.TRUNC){
        //         skipLocus = true;
        //         (*flog) << "Skipping calculation at position " << hm->mapData->mapEntries[locus].physicalPos << " id: " << hm->mapData->mapEntries[locus].locusName;
        //     }
        //     (*flog) << "\n";
        //     pthread_mutex_unlock(&mutex_log);
        //     break;
        // }
        // else if ( physicalDistance(nextLocus, downstream) > p.MAX_GAP )  
        // {
        //     pthread_mutex_lock(&mutex_log);
        //     (*flog) << "WARNING: Reached a gap of " << physicalDistance(nextLocus, downstream) << "bp > " << p.MAX_GAP 
        //     << "bp. Skipping calculation at position " <<  hm->mapData->mapEntries[locus].physicalPos << " id: " <<  hm->mapData->mapEntries[locus].locusName << "\n";
        //     pthread_mutex_unlock(&mutex_log);
        //     skipLocus = true;
        //     break;
        // }
    }
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

    
    double ihh1=0;
    double ihh0=0;

    unsigned int n_c0= hm->hapData->get_n_c0(locus);
    unsigned int n_c1= hm->hapData->get_n_c1(locus);

    uint64_t prev_ehh0_before_norm = 0;
    uint64_t prev_ehh1_before_norm = 0;

    uint64_t curr_ehh0_before_norm = 0;
    uint64_t curr_ehh1_before_norm = 0;

    const double& normalizer_0 = twice_num_pair_or_square(n_c0, p.ALT);
    const double& normalizer_1 = twice_num_pair_or_square(n_c1, p.ALT);

    int numSnps = hm->hapData->nloci;
    int numHaps = hm->hapData->nhaps;

    // PHASE 1: INITIALIZATION
    int* group_count = new int[numHaps];
    int* group_id = new int[numHaps];
    bool* isDerived = new bool [numHaps];
    bool* isAncestral = new bool [numHaps];

    for(int i = 0; i<numHaps; i++){ //will be vectorized with compile time flags
        group_count[i] = 0;
        group_id[i] = 0;
        isDerived[i] = false;
        isAncestral[i] = false;        //assert(hm->hapData->hapEntries[locus].flipped == false);
    }
    int totgc = 0; // total group count upto this point

    // PHASE 1a: INIT CORE LOCUS
    if(n_c1==0){    // all 0s
        group_count[0] = numHaps;
        totgc+=1;
        prev_ehh0_before_norm = normalizer_0;
        bool skipLocus = true;
    }else if (n_c1==numHaps){ // all 1s
        group_count[0] = numHaps;
        totgc+=1;
        if(p.LOW_MEM){
            MyBitset* vb = hm->hapData->hapEntries[locus].hapbitset;
            ACTION_ON_ALL_SET_BITS(vb, {
                isDerived[set_bit_pos] = true;
            });
        }else{
            for (const unsigned int& set_bit_pos : hm->hapData->hapEntries[locus].positions){
                isDerived[set_bit_pos] = true;
            }
        }
        bool skipLocus = true;
        
        prev_ehh1_before_norm = normalizer_1;
    }else{  //so both n_c1 and n_c0 is non-0
        group_count[1] = n_c1;
        group_count[0] = n_c0;

        if(p.LOW_MEM){
            ACTION_ON_ALL_SET_BITS(hm->hapData->hapEntries[locus].hapbitset, {
                isDerived[set_bit_pos] = true;
                group_id[set_bit_pos] = 1;
            });
        }else{
            for (const unsigned int& set_bit_pos : hm->hapData->hapEntries[locus].positions){
                isDerived[set_bit_pos] = true;
                group_id[set_bit_pos] = 1;
            }
        }
        
        totgc+=2;
        prev_ehh0_before_norm = normalizer_0;
        prev_ehh1_before_norm = normalizer_1;
    }

    // PHASE 2: IHS LOOP
    
    if(downstream){
        if(normalizer_1!=0){
            ihh1 += (curr_ehh1_before_norm + prev_ehh1_before_norm) * 0.5 / normalizer_1;
        }
        if(normalizer_0!=0){
            ihh0 += (curr_ehh0_before_norm + prev_ehh0_before_norm) * 0.5 /  normalizer_0;
        }
    }
    
    curr_ehh1_before_norm = prev_ehh1_before_norm;
    curr_ehh0_before_norm = prev_ehh0_before_norm;

    int i = locus;  
    while(true){ // Upstream: for ( int i = locus+1; i<all_positions.size(); i++ )
        if(downstream){
            if (--i < 0) {
                std::cerr<<"Break reason for locus "<<locus<<":: REACHED_LEFT_EDGE."<<endl;
                break;
            }
        }else{
            if (++i >= numSnps) {
                std::cerr<<"Break reason for locus "<<locus<<":: REACHED_RIGHT_EDGE."<<endl;
                break;
            }
        }
        
        if(curr_ehh1_before_norm*1.0/twice_num_pair_or_square(n_c1, p.ALT) <= p.EHH_CUTOFF and curr_ehh0_before_norm*1.0/twice_num_pair_or_square(n_c0, p.ALT)  <= p.EHH_CUTOFF){   // or cutoff, change for benchmarking against hapbin
            //std::cerr<<"Break reason for locus "<<locus<<":: EHH_CUTOFF."<<endl;
            break;
        }

        double distance =  downstream? hm->mapData->mapEntries[i+1].physicalPos - hm->mapData->mapEntries[i].physicalPos : hm->mapData->mapEntries[i].physicalPos - hm->mapData->mapEntries[i-1].physicalPos;
        
        if (distance < 0) // this should not be needed as we already did integrity check previously
        {
            std::cerr << "ERROR: physical position not in ascending order.\n"; 
            throw 0;
        }
        
        if(distance > p.SCALE_PARAMETER){
            distance /= p.SCALE_PARAMETER;
        }
        if(p.CALC_NSL){
            distance = 1;
        }
        
        if(hm->hapData->get_n_c0(i) == 0 or hm->hapData->get_n_c1(i) == 0 ){ // monomorphic check
            pthread_mutex_lock(&mutex_log);
            std::cerr<<"ERROR: Monomorphic site should not exist."<<endl;
            std::cerr<< hm->hapData->get_n_c0(i) <<" n_c0 at locus "<< i <<endl; 
            std::cerr<< hm->hapData->get_n_c1(i) <<" n_c1 at locus "<< i<< endl; 
            pthread_mutex_unlock(&mutex_log);
            throw 0;
            
            if(normalizer_1!=0){    // case where not all are 0s
                ihh1 += (curr_ehh1_before_norm + prev_ehh1_before_norm) * distance * 0.5 / normalizer_1;
            }
            if(normalizer_0!=0){   // case where not all are 1s
                ihh0 += (curr_ehh1_before_norm + prev_ehh1_before_norm) * distance * 0.5 / normalizer_0;
            }
            continue;
        }

        // PHASE 2: GET MAP
        std::unique_ptr<std::unordered_map<unsigned int, std::vector<unsigned int>>> mp(new std::unordered_map<unsigned int, std::vector<unsigned int>>());
        //unordered_map<unsigned int, vector<unsigned int> >  * mp  = new unordered_map<unsigned int, vector<unsigned int> >();
        unordered_map<unsigned int, vector<unsigned int> >& m = (* mp);


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
            vector<unsigned int> * v2p = downstream? &(hm->hapData->hapEntries[i+1].xors) : &(hm->hapData->hapEntries[i].xors);

             // ensure that in boundary we don't do any xor calculation
            if( (hm->hapData->hapEntries[i].positions.size() <(*v2p).size() && i!=numHaps-1) ||  hm->hapData->benchmark_flag != "XOR"){ 
                v2p = &hm->hapData->hapEntries[i].positions;
            }
            //v2p = &hm->hapData->hapEntries[i].positions; // uncomment to disable xor

            for (const unsigned int& set_bit_pos : (*v2p)){
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
            
            bool isDerivedGroup =  (!hm->hapData->hapEntries[locus].flipped && isDerived[ele.second[0]]) || (hm->hapData->hapEntries[locus].flipped && !isAncestral[ele.second[0]]); // just check first element to know if it is derived. 
            if(isDerivedGroup) // if the core locus for this chr has 1 (derived), then update ehh1, otherwise ehh0
            {
                prev_ehh1_before_norm += del_update;
            }else{
                prev_ehh0_before_norm += del_update;
            }
        }

        // CLEAR THE MAP
        m.clear();
        //unordered_map<unsigned int, vector<unsigned int> >().swap(m);

        if(normalizer_1!=0){
            if(prev_ehh1_before_norm*1.0/normalizer_1 > p.EHH_CUTOFF){
                ihh1 += (curr_ehh1_before_norm + prev_ehh1_before_norm) * distance * 0.5 / normalizer_1;
            }
        }

        if(normalizer_0!=0){
            if(prev_ehh0_before_norm*1.0/twice_num_pair_or_square(n_c0, p.ALT) > p.EHH_CUTOFF){
                ihh0 += (curr_ehh0_before_norm + prev_ehh0_before_norm) * distance * 0.5 / normalizer_0;
            }
        }

        curr_ehh1_before_norm = prev_ehh1_before_norm;
        curr_ehh0_before_norm = prev_ehh0_before_norm;

        if(totgc == numHaps) {
            std::cerr<<"Break reason for locus "<<locus<<":: ALL_UNIQUE."<<endl;
            break;
        }
        if(!p.CALC_NSL && physicalDistance_from_core(i,locus, downstream) >= max_extend) {
            std::cerr<<"Break reason for locus "<<locus<<":: MAX_EXTEND."<<endl;
            break;
        }
        if(p.CALC_NSL && abs(i-locus) >= p.MAX_EXTEND_NSL) {
            std::cerr<<"Break reason for locus "<<locus<<":: MAX_EXTEND_NSL."<<endl;
            break; //g(xi−1, xi) = 1.
        }
    }

    delete[] group_count;
    delete[] group_id;
    delete[] isDerived;
    delete[] isAncestral;
    //delete mp;
    return make_pair(ihh1, ihh0);
}


// void IHS::thread_ihs(int tid,  IHS* ehh_obj){
//     // another way to parallelize
//     int thread_id = tid;
//     for (int i = thread_id ; i < ehh_obj->hm->hapData->nloci; i += ehh_obj->numThreads) {
//         ehh_obj->calc_ihh(i);
//     }
//     pthread_mutex_lock(&mutex_log);
//     (*(ehh_obj->flog))<<("Finishing thread # "+to_string(tid)+" at "+to_string(MainTools::readTimer())+"\n");
//     pthread_mutex_unlock(&mutex_log);
// }


void IHS::thread_ihs(int tid,   IHS* ehh_obj){
    int elem_per_block = floor(ehh_obj->hm->hapData->nloci/ehh_obj->numThreads);
    int start = tid*elem_per_block ;
    int end = start + elem_per_block  ;
    if(tid == ehh_obj->numThreads-1 ){
        end = ehh_obj->hm->hapData->nloci;
    }
    //#pragma omp parallel 
    for(int locus = start; locus< end; locus++){
        auto ihh1_ihh0 = ehh_obj->calc_ihh(locus);
        ehh_obj->iHH1[locus] = ihh1_ihh0.first;
        ehh_obj->iHH0[locus] = ihh1_ihh0.second;
    }
    
    pthread_mutex_lock(&mutex_log);
    (*(ehh_obj->flog))<<("Finishing thread # "+to_string(tid)+" at "+to_string(SelscanStats::readTimer())+"\n");
    pthread_mutex_unlock(&mutex_log);
}

// pair<double, double> IHS::thread_ihs(int locus,  unordered_map<unsigned int, vector<unsigned int> >& m, IHS* ehh_obj){
//     // unordered_map<unsigned int, vector<unsigned int> > m2;
//     // m2.clear();
//     return ehh_obj->calc_ihh(locus);
    
// }

void IHS::main() {
    // iHH0 = new double[hm->hapData->nloci];
    // iHH1 = new double[hm->hapData->nloci];
    //std::vector< std::unordered_map<unsigned int, std::vector<unsigned int> > > map_per_thread( hm->hapData->nloci );

    if(hm->hapData->unphased){
        iHH2 = new double[hm->hapData->nloci];
        ciHH0 = new double[hm->hapData->nloci];
        ciHH1 = new double[hm->hapData->nloci];
        ciHH2 = new double[hm->hapData->nloci];
    }
    int total_calc_to_be_done = hm->hapData->nloci;
    cout<<"total calc to be done: "<<total_calc_to_be_done<<endl;

    ThreadPool pool(p.numThreads);
        std::vector< std::future<pair<double, double> > > results;
        for(int i = 0; i <  hm->mapData->nloci; ++i) {
            results.emplace_back(
                pool.enqueue([i,this] {
                    //ihsfinder.calc_ehh_unidirection(i, false);
                    return this->calc_ihh1(i);
                    //return log10(ihsfinder.iHH1[i]/ihsfinder.iHH0[i]);
                    //return 0.0;
                })
            );
            //std::this_thread::sleep_for(std::chrono::milliseconds(1));  // Sleep for 1 second

            
        }
        int locus = 0;
        for(auto && result: results){
            //result.get();
            //std::cout << result.get() << ' ';
            pair<double, double> ihh1_ihh0 = result.get();
            double ihh1 = ihh1_ihh0.first;
            double ihh0 = ihh1_ihh0.second;

            (*fout) << std::fixed <<   hm->mapData->mapEntries[locus].locusName << " " <<   hm->mapData->mapEntries[locus].physicalPos << " "
                        <<  hm->mapData->mapEntries[locus].locId << " " << calcFreq(locus) << " "
                        << ihh1 << " " << ihh0 <<" "<< log10(ihh1/ihh0) <<endl;

            locus++;
            assert(locus<=hm->mapData->nloci);
        }


        /*   
        if(hm->hapData->unphased){
            // for (int locus = 0; locus < hm->hapData->nloci; locus++){
            //     (*fout) << std::fixed <<   hm->mapData->mapEntries[locus].locusName << " " <<   hm->mapData->mapEntries[locus].physicalPos << " "
            //             <<  hm->mapData->mapEntries[locus].locId << " " << hm->hapData->calcFreq(locus) << " "
            //             << iHH2[locus] << " " << iHH0[locus] <<" "<< get_ihs_unphased(locus) <<endl;
            //             //<< log10(ciHH2[locus]/iHH2[locus]) << " " << log10(ciHH0[locus]/iHH0[locus]) <<" "<< get_ihs_unphased(locus) <<endl;
            // }
            // delete[] iHH2;
            // delete[] ciHH0;
            // delete[] ciHH1;
            // delete[] ciHH2;
        }else{
            for (int locus = 0; locus < hm->hapData->nloci; locus++){
                (*fout) << std::fixed <<   hm->mapData->mapEntries[locus].locusName << " " <<   hm->mapData->mapEntries[locus].physicalPos << " "
                        <<  hm->mapData->mapEntries[locus].locId << " " << hm->hapData->calcFreq(locus) << " "
                        << iHH1[locus] << " " << iHH0[locus] <<" "<< log10(iHH1[locus]/iHH0[locus]) <<endl;
                        //<< log10(ciHH2[locus]/iHH2[locus]) << " " << log10(ciHH0[locus]/iHH0[locus]) <<" "<< get_ihs_unphased(locus) <<endl;
            }
        }
        */


        //std::cout << std::endl;
        //print_results_to_out();
    
    //this->runTasks();
    //this->runTaskTQ();
    
   


    //ThreadPool pool(p.numThreads); // Create a thread pool with 4 threads
    
    
    
    
    // for (int i = 0; i < 10; ++i) {
    //     //pool.enqueue([this, i] { calc_ihh(i, map_per_thread[i]) });
    //      pool.enqueue([this, i, map_per_thread[i] ] { calc_ihh(i, map_per_thread[i]); });
    // }}
    
    
    
    // for (int i = 0; i < hm->mapData->nloci; ++i) {
    //   //  pool.enqueue([this, i](std::unordered_map< unsigned int, std::vector<unsigned int> >& map_per_worker) { IHS::thread_ihs(i, map_per_worker, this); });
    //    // pool.enqueue([i](std::unordered_map< unsigned int, std::vector<unsigned int> >& map_per_worker) { thread_ihs(i, map_per_worker, this); });
    //       pool.enqueue(std::bind(&IHS::thread_ihs, i, std::ref(map_per_thread[i]), this));
    // }
    // Give some time for tasks to complete
    //std::this_thread::sleep_for(std::chrono::seconds(2));
    
   hm->hapData->releaseHapData_bitset();
}

void IHS::main_old(){
    iHH0 = new double[hm->hapData->nloci];
    iHH1 = new double[hm->hapData->nloci];

    if(hm->hapData->unphased){
        iHH2 = new double[hm->hapData->nloci];
        ciHH0 = new double[hm->hapData->nloci];
        ciHH1 = new double[hm->hapData->nloci];
        ciHH2 = new double[hm->hapData->nloci];
    }
    
    bool openmp_enabled = false;
    // two different ways to parallelize: first block does pthread, second block does openmp
    if (!openmp_enabled)
    {
        int total_calc_to_be_done = hm->hapData->nloci;
        cout<<"total calc to be done: "<<total_calc_to_be_done<<endl;


        ///*old wrong
        std::thread *myThreads = new std::thread[numThreads];
        for (int i = 0; i < numThreads; i++)
        {
            myThreads[i] = std::thread(thread_ihs, i,  this);
            //myThreads[i] = std::thread(thread_ihs, i, std::ref(map_per_thread[i]), this);
        }
        for (int i = 0; i < numThreads; i++)
        {
            myThreads[i].join(); // Join will block our main thread, and so the program won't exit until all finish
        }
        delete[] myThreads;
        //*/

        (*flog)<<("all threads finished. now calculating ihh...\n");
    }else{
        /*
        // #pragma clang loop unroll_count(8) // 
        // #pragma clang loop vectorize(assume_safety)
        (*flog)<<("open mp enabled. "+to_string(omp_get_max_threads())+" threads\n");
        

        #pragma omp parallel shared(hm)
        {
            #pragma omp for schedule(dynamic,10)
            for(int i = 0 ; i< numSnps; i++){
                calc_ehh(i);
                //cout<<"open mp enabled. "<<omp_get_num_threads()<<"threads"<<endl;
            }
        }
        (*flog)<<("finishing all threads # at "+to_string(readTimer())+"\n");
        */
    }
   

    /*
    char str[80];
    for (int i = 0; i < numSnps; i++){
         if(hm.getMAF(i) >= min_maf && 1-hm.getMAF(i) >= min_maf){
            sprintf(str, "%d %d %f %f %f %f\n", hm.mentries[i].phyPos, hm.mentries[i].locId, hm.all_positions[i].size()*1.0/numHaps, iHH1[i], iHH0[i], log10(iHH1[i]/ iHH0[i]));
            out_ihs<<str;
         }
    }
    */

   
    if(hm->hapData->unphased){
        for (int locus = 0; locus < hm->hapData->nloci; locus++){
            (*fout) << std::fixed <<   hm->mapData->mapEntries[locus].locusName << " " <<   hm->mapData->mapEntries[locus].physicalPos << " "
                    <<  hm->mapData->mapEntries[locus].locId << " " << calcFreq(locus) << " "
                    << iHH2[locus] << " " << iHH0[locus] <<" "<< get_ihs_unphased(locus) <<endl;
                    //<< log10(ciHH2[locus]/iHH2[locus]) << " " << log10(ciHH0[locus]/iHH0[locus]) <<" "<< get_ihs_unphased(locus) <<endl;
        }
        delete[] iHH2;
        delete[] ciHH0;
        delete[] ciHH1;
        delete[] ciHH2;
    }else{
        for (int locus = 0; locus < hm->hapData->nloci; locus++){
            (*fout) << std::fixed <<   hm->mapData->mapEntries[locus].locusName << " " <<   hm->mapData->mapEntries[locus].physicalPos << " "
                    <<  hm->mapData->mapEntries[locus].locId << " " << calcFreq(locus) << " "
                    << iHH1[locus] << " " << iHH0[locus] <<" "<< log10(iHH1[locus]/iHH0[locus]) <<endl;
                    //<< log10(ciHH2[locus]/iHH2[locus]) << " " << log10(ciHH0[locus]/iHH0[locus]) <<" "<< get_ihs_unphased(locus) <<endl;
        }
    }
    delete[] iHH0;
    delete[] iHH1; 
    hm->hapData->releaseHapData_bitset();

    // if(hm->hapData->unphased){
    //     (*fout) << std::fixed <<   hm->mapData->mapEntries[locus].locusName << " " <<   hm->mapData->mapEntries[locus].physicalPos << " "
    //             <<  hm->mapData->mapEntries[locus].locId << " " << hm->hapData->calcFreq(locus) << " "
    //             << iHH1[locus] << " " << iHH0[locus] <<" "<< get_ihs_unphased(locus) <<endl;
    // }

    return;
}


/**
 * @brief Calculate unphased IHS
 * @param locus The locus 
*/
double IHS::get_ihs_unphased(int locus){
                //  ihh1[locus] = log10(derived_ihh / notDerived_ihh);
                // ihh2[locus] = log10(ancestral_ihh / notAncestral_ihh);
                // ihs[locus] = (ihh1[locus] > ihh2[locus]) ? ihh1[locus] : 0-ihh2[locus];
     double iHS2 = log10(iHH2[locus]/ciHH2[locus]);
     double iHS0 = log10(iHH0[locus]/ciHH0[locus]);
    if(iHS2 > iHS0){
        return iHS2;
    }else{
        return 0-iHS0;
    }
}



/**
 * @brief Calculate the EHH for a single locus as part of IHH routine
 * @param locus The locus 
*/
pair<double, double> IHS::calc_ihh1(int locus){
    // unordered_map<unsigned int, vector<unsigned int> >* mp = new unordered_map<unsigned int, vector<unsigned int> >();
    // unordered_map<unsigned int, vector<unsigned int> > m = *mp;
    int numSnps = hm->mapData->nloci;
    int numHaps = hm->hapData->nhaps;
    
    //unordered_map<unsigned int, vector<unsigned int> > m;
    // iHH0[locus] = 0;
    // iHH1[locus] = 0;

    // double& ihh1=iHH1[locus];
    // double& ihh0=iHH0[locus];

    double ihh1;
    double ihh0;

    if(hm->hapData->unphased){
        iHH2[locus] = 0;
        ciHH0[locus] = 0;
        ciHH1[locus] = 0;
        ciHH2[locus] = 0;
        calc_ehh_unidirection_unphased(locus, false); // upstream
        calc_ehh_unidirection_unphased(locus, true); // downstream
    }else{

        auto ihh1_ihh0_upstream = calc_ehh_unidirection(locus, false); // upstream
            auto ihh1_ihh0_downstream = calc_ehh_unidirection(locus,  true); // downstream
            ihh1 = ihh1_ihh0_upstream.first + ihh1_ihh0_downstream.first;
            ihh0 = ihh1_ihh0_upstream.second + ihh1_ihh0_downstream.second;
        // if(p.LOW_MEM){
        //     auto ihh1_ihh0_upstream = calc_ehh_unidirection_bitset(locus, false); // upstream
        //     auto ihh1_ihh0_downstream = calc_ehh_unidirection_bitset(locus, true); // downstream
        //     ihh1 = ihh1_ihh0_upstream.first + ihh1_ihh0_downstream.first;
        //     ihh0 = ihh1_ihh0_upstream.second + ihh1_ihh0_downstream.second;
        // }else{
        //     auto ihh1_ihh0_upstream = calc_ehh_unidirection(locus, false); // upstream
        //     auto ihh1_ihh0_downstream = calc_ehh_unidirection(locus,  true); // downstream
        //     ihh1 = ihh1_ihh0_upstream.first + ihh1_ihh0_downstream.first;
        //     ihh0 = ihh1_ihh0_upstream.second + ihh1_ihh0_downstream.second;
        // }
    }
    
    if(!hm->hapData->unphased){
        //handle all corner cases
        if(!p.LOW_MEM){
            if(hm->hapData->hapEntries[locus].positions.size()==0){
                ihh1 = 1;
            }
            if(hm->hapData->hapEntries[locus].positions.size()==numHaps){ 
                ihh0 = 1;
            }
            if(hm->hapData->hapEntries[locus].positions.size()==1){
                ihh1 = (locus == 0 or locus == numSnps-1)? 0.5 : 1;
            }
            if(hm->hapData->hapEntries[locus].positions.size()==numHaps-1){
                ihh0 = (locus == 0 or locus == numSnps-1)? 0.5 : 1;
            }
            //PRINTHERE
            //  if(hm.getMAF(i) >= min_maf && 1-hm.getMAF(i) >= min_maf){
            //         sprintf(str, "%d %d %f %f %f %f\n", hm.mentries[i].phyPos, hm.mentries[i].locId, hm.all_positions[i].size()*1.0/numHaps, iHH1[i], iHH0[i], log10(iHH1[i]/ iHH0[i]));
            //         out_ihs<<str;
            //      }

            // (*fout) << std::fixed <<   hm->mapData->mapEntries[locus].locusName << " " <<   hm->mapData->mapEntries[locus].physicalPos << " "
            //         <<  hm->mapData->mapEntries[locus].locId << " " << hm->hapData->calcFreq(locus) << " "
            //         << iHH1[locus] << " " << iHH0[locus] <<" "<< iHH1[locus]*1.0/iHH0[locus] <<endl;
        
        }
    }
    //delete mp;
    return make_pair(ihh1, ihh0);
}



/**
 * @brief Calculate the EHH for a single locus as part of IHH routine
 * @param locus The locus 
*/
pair<double, double>  IHS::calc_ihh(int locus){
    // unordered_map<unsigned int, vector<unsigned int> >* mp = new unordered_map<unsigned int, vector<unsigned int> >();
    // unordered_map<unsigned int, vector<unsigned int> > m = *mp;
    int numSnps = hm->mapData->nloci;
    int numHaps = hm->hapData->nhaps;
    
    //unordered_map<unsigned int, vector<unsigned int> > m;
    // iHH0[locus] = 0;
    // iHH1[locus] = 0;

    // double& ihh1=iHH1[locus];
    // double& ihh0=iHH0[locus];
    double ihh1, ihh0;

    if(hm->hapData->unphased){
        iHH2[locus] = 0;
        ciHH0[locus] = 0;
        ciHH1[locus] = 0;
        ciHH2[locus] = 0;
        calc_ehh_unidirection_unphased(locus, false); // upstream
        calc_ehh_unidirection_unphased(locus, true); // downstream
    }else{
        if(p.LOW_MEM){
            //calc_ehh_unidirection_bitset(locus, false, m); // upstream
            //calc_ehh_unidirection_bitset(locus, true, m); // downstream
        }else{
            // calc_ehh_unidirection(locus, false); // upstream
            // calc_ehh_unidirection(locus,  true); // downstream
            auto ihh1_ihh0_upstream = calc_ehh_unidirection(locus, false); // upstream
            auto ihh1_ihh0_downstream = calc_ehh_unidirection(locus,  true); // downstream
            ihh1 = ihh1_ihh0_upstream.first + ihh1_ihh0_downstream.first;
            ihh0 = ihh1_ihh0_upstream.second + ihh1_ihh0_downstream.second;
        }
    }
    
    if(!hm->hapData->unphased){
        //handle all corner cases
        if(hm->hapData->benchmark_flag!="BITSET"){
            if(hm->hapData->hapEntries[locus].positions.size()==0){
                ihh1 = 1;
            }
            if(hm->hapData->hapEntries[locus].positions.size()==numHaps){ 
                ihh0 = 1;
            }
            if(hm->hapData->hapEntries[locus].positions.size()==1){
                ihh1 = (locus == 0 or locus == numSnps-1)? 0.5 : 1;
            }
            if(hm->hapData->hapEntries[locus].positions.size()==numHaps-1){
                ihh0 = (locus == 0 or locus == numSnps-1)? 0.5 : 1;
            }
            //PRINTHERE
            //  if(hm.getMAF(i) >= min_maf && 1-hm.getMAF(i) >= min_maf){
            //         sprintf(str, "%d %d %f %f %f %f\n", hm.mentries[i].phyPos, hm.mentries[i].locId, hm.all_positions[i].size()*1.0/numHaps, iHH1[i], iHH0[i], log10(iHH1[i]/ iHH0[i]));
            //         out_ihs<<str;
            //      }

            // (*fout) << std::fixed <<   hm->mapData->mapEntries[locus].locusName << " " <<   hm->mapData->mapEntries[locus].physicalPos << " "
            //         <<  hm->mapData->mapEntries[locus].locId << " " << hm->hapData->calcFreq(locus) << " "
            //         << iHH1[locus] << " " << iHH0[locus] <<" "<< iHH1[locus]*1.0/iHH0[locus] <<endl;
        
        }
    }
    //delete mp;
    return make_pair(ihh1, ihh0);
}


