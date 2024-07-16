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
        
        //double del_update = -twice_num_pair_or_square(group_count[old_group_id], p.ALT) + twice_num_pair_or_square(newgroup_size, p.ALT) + twice_num_pair_or_square(group_count[old_group_id] - newgroup_size, p.ALT);
        
        group_count[old_group_id] -= newgroup_size;
        group_count[totgc] += newgroup_size;
        
        totgc+=1;
        
        //TODO
        if(is1[ele.second[0]]) // if the core locus for this chr has 1 (derived), then update ehh1, otherwise ehh0
        {
            group_core[totgc-1] = 1;
            //ehh_before_norm[1] += del_update;
            
            // cehh_before_norm[2] += del_update;
            // cehh_before_norm[0] += del_update;
        }else if(is2[ele.second[0]]) {
            group_core[totgc-1] = 2;
            //ehh_before_norm[2] += del_update;

            // cehh_before_norm[1] += del_update;
            // cehh_before_norm[0] += del_update;
        }else{
            group_core[totgc-1] = 0;

            //ehh_before_norm[0] += del_update;

            // cehh_before_norm[2] += del_update;
            // cehh_before_norm[1] += del_update;
        }

        
    }
}

pair<double, double> IHS::calc_ehh_unidirection_unphased(int locus, bool downstream,  double& cihh2, double& cihh0){
    bool DEBUG_UNP = true;

    std::unique_ptr<std::unordered_map<int, std::vector<int>>> mp(new std::unordered_map<int, std::vector<int>>());
    unordered_map<int, vector<int> >& m2 = (* mp);

    int numSnps = hm->hapData->nloci;
    int numHaps = hm->hapData->nhaps;

    double prev_ehh_before_norm[3] = {0,0,0};
    double curr_ehh_before_norm[3] = {0,0,0};

    double  prev_cehh_before_norm[3] = {1,1,1};
    double  curr_cehh_before_norm[3] = {1,1,1};

    int n_c[3] = {0,0,0};

    int* group_count = new int[numHaps];
    int* group_id =  new int[numHaps];
    int* group_core =  new int[numHaps];

    bool* is1 =  new bool[numHaps];
    bool* is2 =  new bool[numHaps];

//     std::unique_ptr<int[]> group_count = std::make_unique<int[]>(numHaps);
// std::unique_ptr<int[]> group_id = std::make_unique<int[]>(numHaps);
// std::unique_ptr<bool[]> is1 = std::make_unique<bool[]>(numHaps);
// std::unique_ptr<bool[]> is2 = std::make_unique<bool[]>(numHaps);

    double iHH[3] = {0,0,0};
    double ciHH[3] = {0,0,0};
    // double* iHH[3] = {iHH0, iHH1, iHH2};
    // double* ciHH[3] = {ciHH0, ciHH1, ciHH2};

    //will be vectorized with compile time flags
    for(int i = 0; i<numHaps; i++){
        group_count[i] = 0;
        is1[i] = false;
        is2[i] = false;
    }


    int totgc=0;

    if(p.LOW_MEM){
        // n_c[0] = numHaps - hm->hapData->hapEntries[locus].xorbitset->num_1s;
        // n_c[1] = hm->hapData->hapEntries[locus].hapbitset->num_1s;
        // n_c[2] = numHaps-n_c[0]-n_c[1];
        // // // DEBUG cout<<"n_c0: "<<n_c[0]<<" n_c1: "<<n_c[1]<<" n_c2: "<<n_c[2]<<endl;

       
        n_c[1] = hm->hapData->hapEntries[locus].hapbitset->num_1s;
        n_c[2] = hm->hapData->hapEntries[locus].xorbitset->num_1s;
        n_c[0] = numHaps - n_c[1] - n_c[2];
    }else{
        n_c[1] = hm->hapData->hapEntries[locus].positions.size();
        n_c[2] = hm->hapData->hapEntries[locus].positions2.size();
        n_c[0] = numHaps - n_c[1] - n_c[2];
        // // // DEBUG cout<<"n_c0: "<<n_c[0]<<" n_c1: "<<n_c[1]<<" n_c2: "<<n_c[2]<<endl;
    }


    if(n_c[1] + n_c[2] + n_c[0] != numHaps){
        cerr<<"ERROR: n_c1 + n_c2 + n_c0 != numHaps"<<endl;
        exit(1);
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

    // if(locus==0){
    //     cout<<"normalizer: "<<normalizer[0]<<" "<<normalizer[1]<<" "<<normalizer[2]<<endl;
    //     cout<<"normalizer_not: "<<normalizer_not[0]<<" "<<normalizer_not[1]<<" "<<normalizer_not[2]<<endl;
    // }    

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

   if(p.LOW_MEM){
        ACTION_ON_ALL_SET_BITS(hm->hapData->hapEntries[locus].hapbitset, {
            is1[set_bit_pos] = true;
            group_id[set_bit_pos] = pos_of_012[1];   
        });

        ACTION_ON_ALL_SET_BITS(hm->hapData->hapEntries[locus].xorbitset, {
            is2[set_bit_pos] = true;
            group_id[set_bit_pos] = pos_of_012[2]; 
        });
                    // if(is1[set_bit_pos] == false){
            //     is2[set_bit_pos] = true;
            //     group_id[set_bit_pos] = pos_of_012[2];         
            // }    

   }else{
        const vector<unsigned int>& v2 = hm->hapData->hapEntries[locus].positions2;
        for(const int& set_bit_pos : v2){
            is2[set_bit_pos] = true;
            group_id[set_bit_pos] = pos_of_012[2];
        }

        const vector<unsigned int>& v1 = hm->hapData->hapEntries[locus].positions;
        for(const int& set_bit_pos : v1){
            is1[set_bit_pos] = true;
            group_id[set_bit_pos] = pos_of_012[1];
        }
   }

    // // // DEBUG::: cout<<locus<<" :: ";
    // // // DEBUG:::for (int i = 0; i<numHaps; i++){
    // // // DEBUG:::    if(group_id[i]!=0 and group_id[i]!=1)
    // // // DEBUG:::         cout<<i<<"->"<<group_id[i]<<". ";
    // // // DEBUG:::}
    // // // DEBUG::: cout<<endl;

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

    // if(downstream){
    //     for(int i=0; i<3; i++){
    //         //we dont need i = 1
    //         if(twice_num_pair_or_square(n_c[i], p.ALT)!=0){
    //             iHH[i] += (curr_ehh_before_norm[i] + prev_ehh_before_norm[i]) * 0.5 / twice_num_pair_or_square(n_c[i], p.ALT);
    //         }

    //         if(twice_num_pair_or_square(numHaps - n_c[i], p.ALT)!=0){
    //             ciHH[i] += (curr_cehh_before_norm[i] + prev_cehh_before_norm[i])  * 0.5 / twice_num_pair_or_square(numHaps - n_c[i], p.ALT);   
    //         }
    //     }
    // }
    
    for(int i=0; i<3; i++){
        //we dont need i = 1
        prev_ehh_before_norm[i] = curr_ehh_before_norm[i];
        prev_cehh_before_norm[i] =  curr_cehh_before_norm[i];
    }

    // if(DEBUG_UNP){
    //     if(hm->mapData->mapEntries[locus].physicalPos == 2288883){
    //         cout<<"locus 1 "<<locus<<" "<<hm->hapData->hapEntries[locus].positions.size()<<endl;
    //         cout<<"locus 2 "<<locus<<" "<<hm->hapData->hapEntries[locus].positions2.size()<<endl;
    //         cout<<"locus 0 "<<locus<<" "<<numHaps- hm->hapData->hapEntries[locus].positions.size() -hm->hapData->hapEntries[locus].positions2.size()<<endl;
    //         exit(1);
    //     }
    // }

    int i = locus;  
    while(true){ // Upstream: for ( int i = locus+1; i<all_positions.size(); i++ )
        if(curr_ehh_before_norm[2]*1.0/normalizer[2] <= p.EHH_CUTOFF and curr_ehh_before_norm[0]*1.0/normalizer[0]  <= p.EHH_CUTOFF){   // or cutoff, change for benchmarking against hapbin
            //std::cerr<<"Break reason for locus "<<locus<<":: EHH_CUTOFF."<<endl;
            break;
        }

        bool edgeBreak = false;
        edgeBreak = (downstream)? (i-1 < 0) : (i+1 >= numSnps);
        if(edgeBreak) {
            pthread_mutex_lock(&mutex_log);
            (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << p.EHH_CUTOFF
                    << ". ";
            if (!p.TRUNC){
                bool skipLocus = 1;
                (*flog) << "Skipping calculation at position " << hm->mapData->mapEntries[locus].physicalPos << " id: " << hm->mapData->mapEntries[locus].locusName;
            }
            (*flog) << endl;
            pthread_mutex_unlock(&mutex_log);
        }

        if(downstream){
            --i;
            if(i<0){
                break;
            }
        }else{
            ++i;
            if(i>=numSnps){
                break;
            }
        } 
        
        //long double distance =  downstream? hm->mapData->mapEntries[i+1].geneticPos - hm->mapData->mapEntries[i].geneticPos : hm->mapData->mapEntries[i].geneticPos - hm->mapData->mapEntries[i-1].geneticPos;
        long double distance =  geneticDistance(i, downstream);
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

        if (physicalDistance(i,downstream) > p.MAX_GAP)
        {
            pthread_mutex_lock(&mutex_log);
            (*flog) << "WARNING: Reached a gap of " << physicalDistance(i,downstream)
                    << "bp > " << p.MAX_GAP << "bp. Skipping calculation at position " << hm->mapData->mapEntries[locus].physicalPos << " id: " << hm->mapData->mapEntries[locus].locusName << "\n";
            pthread_mutex_unlock(&mutex_log);
            bool skipLocus = 1;
            break;
        }

        
        if(n_c[0] == numHaps or n_c[1] == numHaps or n_c[2] == numHaps){
            cerr<<"WARNING: Monomorphic site."<<endl;
        
            for(int i=0; i<3; i++){
                if(normalizer[i]!=0){
                    iHH[i] += (curr_ehh_before_norm[i] * 1.0 / normalizer[i]  + prev_ehh_before_norm[i] * 1.0  / normalizer[i] ) * 0.5 * distance;
                }
                if(normalizer_not[i]!=0){
                    ciHH[i] += (curr_cehh_before_norm[i] * 1.0  / normalizer_not[i] + prev_cehh_before_norm[i] * 1.0  / normalizer_not[i]) * 0.5  *  distance;
                }
            }
            continue;
        }

        if(p.LOW_MEM){
            
            //unordered_map<int, vector<int> > m = (* mp);
            //unordered_map<int, vector<int> > m;
            
            // MyBitset* both12_p = hm->hapData->hapEntries[i].xorbitset;
            // MyBitset* only1_p = hm->hapData->hapEntries[i].hapbitset;

            ACTION_ON_ALL_SET_BITS(hm->hapData->hapEntries[i].xorbitset, {
                int old_group_id = group_id[set_bit_pos];
                m2[old_group_id].push_back(set_bit_pos);
            });
            updateEHH_from_split_unphased(m2, group_count, group_id, totgc, curr_ehh_before_norm, curr_cehh_before_norm, is1, is2, group_core);
            m2.clear();

            ACTION_ON_ALL_SET_BITS(hm->hapData->hapEntries[i].hapbitset, {
                int old_group_id = group_id[set_bit_pos];
                m2[old_group_id].push_back(set_bit_pos);
            });
            updateEHH_from_split_unphased(m2, group_count, group_id, totgc, curr_ehh_before_norm, curr_cehh_before_norm, is1, is2, group_core);
            m2.clear();

        }else{
            const vector<unsigned int>& g0 = (downstream)? hm->hapData->hapEntries[i+1].g[0] : hm->hapData->hapEntries[i].g[0];
            const vector<unsigned int>& g1 = (downstream)? hm->hapData->hapEntries[i+1].g[1] : hm->hapData->hapEntries[i].g[1];
            const vector<unsigned int>& g2 = (downstream)? hm->hapData->hapEntries[i+1].g[2] : hm->hapData->hapEntries[i].g[2];

            unordered_map<int, vector<int> > m;
            for (const unsigned int& set_bit_pos : g0){
                int old_group_id = group_id[set_bit_pos];
                m[old_group_id].push_back(set_bit_pos);         
            }
            updateEHH_from_split_unphased(m, group_count, group_id, totgc, curr_ehh_before_norm, curr_cehh_before_norm, is1, is2, group_core);
            m.clear();

            for (const unsigned int& set_bit_pos : g1){
                int old_group_id = group_id[set_bit_pos];
                m[old_group_id].push_back(set_bit_pos);         
            }
            updateEHH_from_split_unphased(m, group_count, group_id, totgc, curr_ehh_before_norm, curr_cehh_before_norm, is1, is2, group_core);
            m.clear();


            for (const unsigned int& set_bit_pos : g2){
                int old_group_id = group_id[set_bit_pos];
                m[old_group_id].push_back(set_bit_pos);         
            }
            updateEHH_from_split_unphased(m, group_count, group_id, totgc, curr_ehh_before_norm, curr_cehh_before_norm, is1, is2, group_core);
            m.clear();
            
            

            
            
 
        }
        DEBUG_UNP = false;
        if(DEBUG_UNP){
            if(i==4){
                cout<<"not 0, not 2 "<< n_c[1]+ n_c[2]<<" "<<n_c[1]+n_c[0]<<endl;
                cout<<"cehh 0 "<< (curr_cehh_before_norm[0] - prev_cehh_before_norm[0])/2<<" cehh 2 "<<(curr_cehh_before_norm[2] - prev_cehh_before_norm[2])/2<<endl;
                for(int j = 0; j<totgc; j++){
                    cout<<j<<" "<< group_count[j] << " " <<group_core[j]<<endl;
                }

                // for(int j = 0; j<numHaps; j++){
                //     if(is2[j]){
                //     cout<<"2 start: " << j<<" "<<group_id[j]<<" "<< group_count[group_id[j]]<<" "<<is1[j]<<" "<<is2[j]<<endl;

                //     }
                // }
                // cout<<"---"<<endl;
                // for(int j = 0; j<numHaps; j++){
                //     if(!is2[j]&& !is1[j]){
                //     cout<<"0 start: " << j<<" "<<group_id[j]<<" "<< group_count[group_id[j]]<<" "<<is1[j]<<" "<<is2[j]<<endl;

                //     }
                // }
            }
            //exit(1);
        }
        double homo2 = 0;
        double homo0 = 0;
        double do2 = 0;
        double do0 = 0;
        for(int x = 0; x<totgc; x++){
            int gcsquare = twice_num_pair_or_square(group_count[x],p.ALT);

            if(group_core[x]==0){
                do0 += gcsquare;
                homo2 += gcsquare;
            }else if(group_core[x]==1){
                homo2 += gcsquare;
                homo0 += gcsquare;
            }else{
                do2 += gcsquare;
                homo0 += gcsquare;
            }
        }

        curr_ehh_before_norm[2] = do2;
        curr_ehh_before_norm[0] = do0;
        curr_cehh_before_norm[2] = homo2;
        curr_cehh_before_norm[0] = homo0;
        //cout<<"cehh 0 "<< (curr_cehh_before_norm[0])<<" cehh 2 "<<(curr_cehh_before_norm[2])<<endl;
        
        

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
           
            if(curr_ehh_before_norm[i]*1.0/normalizer[i] > p.EHH_CUTOFF){
                // this was working iHH[i][locus] += (curr_ehh_before_norm[i] + prev_ehh_before_norm[i]) * scale * distance* 0.5 / twice_num_pair(n_c[i]);
                if(twice_num_pair_or_square(n_c[i], p.ALT)!=0){
                    iHH[i] += (curr_ehh_before_norm[i] * 1.0 / normalizer[i] + prev_ehh_before_norm[i]* 1.0/ normalizer[i])  * 0.5 *distance;
                }


            }
            
            if(normalizer_not[i]!=0){
                //this is how it's in selscan //if(cehh_before_norm[i]*1.0/twice_num_pair(numHaps - n_c[i]) > p.EHH_CUTOFF){
                    ciHH[i] += (curr_cehh_before_norm[i] *0.5 / (normalizer_not[i]) + prev_cehh_before_norm[i] *0.5/ normalizer_not[i])*distance ;
                //}
            }
            prev_ehh_before_norm[i] = curr_ehh_before_norm[i];
            prev_cehh_before_norm[i] = curr_cehh_before_norm[i];
        }


        // if(DEBUG_UNP)
        //     cout<<i<<" :: "<<curr_cehh_before_norm[0]/ (normalizer_not[0]) <<" "<<curr_cehh_before_norm[2]/ (normalizer_not[2])<<endl;
        
                //TODO
            //We've now counted all of the unique haplotypes extending out of the core SNP
            //If locus is monomorphic, shoot a warning and skip locus
            //This probably isn't necessary any more
            // if ( !unphased && (numDerived == 0 || numAncestral == 0) ) 
            // {
            //     pthread_mutex_lock(&mutex_log);
            //     (*flog) << "WARNING: locus " << locusName[locus]
            //             << " (number " << locus + 1 << ") is monomorphic. Skipping calculation at this locus.\n";
            //     pthread_mutex_unlock(&mutex_log);
            //     skipLocus = 1;
            //     break;
            // }


        // double freqHetGT = double(numHet) / double(numDerived + numAncestral + numHet);
        //     //double freqAncestralGT = double(numAncestral) / double(numDerived + numAncestral + numHet);

        //     if ( unphased && freqHetGT > 1-MAF ) 
        //     {
        //         pthread_mutex_lock(&mutex_log);
        //         (*flog) << "WARNING: locus " << locusName[locus]
        //                 << " (number " << locus + 1 << ") has too many hets. Skipping calculation at this locus. "
        //                 << "het: " << numHet << " hom0: " << numAncestral << " hom1: " << numDerived << ".\n";
        //         pthread_mutex_unlock(&mutex_log);
        //         skipLocus = 1;
        //         break;
        //     }

        if (physicalDistance_from_core(i, locus, downstream) >= max_extend) break;

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
    delete[] group_count;
    delete[] group_id;
    delete[] is1;
    delete[] is2;
    delete[] group_core;

    cihh2 = ciHH[2];
    cihh0 = ciHH[0];

    // if(cihh2==0) cihh2=1;
    // if(cihh0==0) cihh0=1;
    // if(iHH[2]==0) iHH[2]=1;
    // if(iHH[0]==0) iHH[0]=1;
    return make_pair(iHH[2], iHH[0]);
    //return make_pair(ciHH[2], ciHH[0]);
    //return make_pair(iHH[2], iHH[0]);


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
                //std::cerr<<"Break reason for locus "<<locus<<":: REACHED_LEFT_EDGE."<<endl;
                break;
            }
        }else{
            if (++i >= numSnps) {
                //std::cerr<<"Break reason for locus "<<locus<<":: REACHED_RIGHT_EDGE."<<endl;
                break;
            }
        }
        
        if(curr_ehh1_before_norm*1.0/twice_num_pair_or_square(n_c1, p.ALT) <= p.EHH_CUTOFF and curr_ehh0_before_norm*1.0/twice_num_pair_or_square(n_c0, p.ALT)  <= p.EHH_CUTOFF){   // or cutoff, change for benchmarking against hapbin
            //std::cerr<<"Break reason for locus "<<locus<<":: EHH_CUTOFF."<<endl;
            break;
        }

        double distance =  physicalDistance(i, downstream);
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
            //std::cerr<<"Break reason for locus "<<locus<<":: ALL_UNIQUE."<<endl;
            break;
        }
        if(!p.CALC_NSL && physicalDistance_from_core(i,locus, downstream) >= max_extend) {
            //std::cerr<<"Break reason for locus "<<locus<<":: MAX_EXTEND."<<endl;
            break;
        }
        if(p.CALC_NSL && abs(i-locus) >= p.MAX_EXTEND_NSL) {
            //std::cerr<<"Break reason for locus "<<locus<<":: MAX_EXTEND_NSL."<<endl;
            break; 
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


void IHS::thread_ihs(int tid,   IHS* ehh_obj, double& iHH1, double& iHH0){
    int elem_per_block = floor(ehh_obj->hm->hapData->nloci/ehh_obj->numThreads);
    int start = tid*elem_per_block ;
    int end = start + elem_per_block  ;
    if(tid == ehh_obj->numThreads-1 ){
        end = ehh_obj->hm->hapData->nloci;
    }
    //#pragma omp parallel 
    for(int locus = start; locus< end; locus++){
        auto ihh1_ihh0 = ehh_obj->calc_ihh1(locus);
        iHH1 = ihh1_ihh0.first;
        iHH0 = ihh1_ihh0.second;
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

    // if(hm->hapData->unphased){
    //     iHH0 = new double[hm->hapData->nloci];
    //     iHH1 = new double[hm->hapData->nloci];
    //     iHH2 = new double[hm->hapData->nloci];
    //     ciHH0 = new double[hm->hapData->nloci];
    //     ciHH1 = new double[hm->hapData->nloci];
    //     ciHH2 = new double[hm->hapData->nloci];
    // }


    int total_calc_to_be_done = hm->hapData->nloci;
    cout<<"Total number of loci: "<<total_calc_to_be_done<<endl;



    if(hm->p.numThreads == 1){
        // for(int locus = 0; locus <  hm->mapData->nloci; locus++) {
        for(int locus = 0; locus <  1; locus++) {
       
            pair<double, double> ihh1_ihh0 = calc_ihh1(locus);
            double ihh1 = ihh1_ihh0.first;
            double ihh0 = ihh1_ihh0.second;

            if(!hm->hapData->unphased){
                (*fout) << std::fixed <<   hm->mapData->mapEntries[locus].locusName << " " <<   hm->mapData->mapEntries[locus].physicalPos << " "
                        <<  hm->mapData->mapEntries[locus].locId << " " << hm->hapData->calcFreq(locus) << " "
                        << ihh1 << " " << ihh0 <<" "<< log10(ihh1/ihh0) <<endl;
            }
            if(hm->hapData->unphased){
                const double& iHS2 = ihh1;
                const double& iHS0 = ihh0;
                double ihs = (iHS2) > (iHS0) ? iHS2 : 0-iHS0;
                (*fout) << std::fixed <<   hm->mapData->mapEntries[locus].locusName << " " <<   hm->mapData->mapEntries[locus].physicalPos << " "
                        <<  hm->mapData->mapEntries[locus].locId << " " << hm->hapData->calcFreq(locus) << " "
                        << iHS2 << " " << iHS0 <<" "<< ihs <<endl;
            }
        }
        
        return;
    }else{
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

            if(!hm->hapData->unphased){
                (*fout) << std::fixed <<   hm->mapData->mapEntries[locus].locusName << " " <<   hm->mapData->mapEntries[locus].physicalPos << " "
                        <<  hm->mapData->mapEntries[locus].locId << " " << hm->hapData->calcFreq(locus) << " "
                        << ihh1 << " " << ihh0 <<" "<< log10(ihh1/ihh0) <<endl;
            }
            //calcFreq(locus)
            
            if(hm->hapData->unphased){
                const double& iHS2 = ihh1;
                const double& iHS0 = ihh0;
                double ihs = iHS2 > iHS0 ? iHS2 : 0-iHS0;
                (*fout) << std::fixed <<   hm->mapData->mapEntries[locus].locusName << " " <<   hm->mapData->mapEntries[locus].physicalPos << " "
                        <<  hm->mapData->mapEntries[locus].locId << " " << hm->hapData->calcFreq(locus) << " "
                        << iHS2 << " " << iHS0 <<" "<< ihs <<endl;
                        //<< log10(ciHH2[locus]/iHH2[locus]) << " " << log10(ciHH0[locus]/iHH0[locus]) <<" "<< get_ihs_unphased(locus) <<endl;
            }
            // delete[] iHH2;
            // delete[] ciHH0;
            // delete[] ciHH1;
            // delete[] ciHH2;

            locus++;
            assert(locus<=hm->mapData->nloci);
        }
    }
   //hm->hapData->releaseHapData_bitset();
}

void IHS::main_old(){

    double ihh1 = 0;
    double ihh0 = 0;

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
            myThreads[i] = std::thread(thread_ihs, i,  this, std::ref(ihh1), std::ref(ihh0));
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
        const double& iHS2 = ihh1;
        const double& iHS0 = ihh0;
        double ihs = iHS2 > iHS0 ? iHS2 : 0-iHS0;
        for (int locus = 0; locus < hm->hapData->nloci; locus++){
            (*fout) << std::fixed <<   hm->mapData->mapEntries[locus].locusName << " " <<   hm->mapData->mapEntries[locus].physicalPos << " "
                    <<  hm->mapData->mapEntries[locus].locId << " " << hm->hapData->calcFreq(locus) << " "
                    << iHS2 << " " << iHS0 <<" "<< ihs <<endl;
                    //<< log10(ciHH2[locus]/iHH2[locus]) << " " << log10(ciHH0[locus]/iHH0[locus]) <<" "<< get_ihs_unphased(locus) <<endl;
        }
    }else{
        for (int locus = 0; locus < hm->hapData->nloci; locus++){
            (*fout) << std::fixed <<   hm->mapData->mapEntries[locus].locusName << " " <<   hm->mapData->mapEntries[locus].physicalPos << " "
                    <<  hm->mapData->mapEntries[locus].locId << " " <<  hm->hapData->calcFreq(locus) << " "
                    << ihh1 << " " << ihh0 <<" "<< log10(ihh1/ihh0) <<endl;
                    //<< log10(ciHH2[locus]/iHH2[locus]) << " " << log10(ciHH0[locus]/iHH0[locus]) <<" "<< get_ihs_unphased(locus) <<endl;
        }
    }

    // if(p.LOW_MEM){
    //     hm->hapData->releaseHapData_bitset();
    // }

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
    cerr<<"get_ihs_unphased not implemented yet"<<endl;
    throw 1;
    //             //  ihh1[locus] = log10(derived_ihh / notDerived_ihh);
    //             // ihh2[locus] = log10(ancestral_ihh / notAncestral_ihh);
    //             // ihs[locus] = (ihh1[locus] > ihh2[locus]) ? ihh1[locus] : 0-ihh2[locus];
    //  double iHS2 = log10(iHH2[locus]/ciHH2[locus]);
    //  double iHS0 = log10(iHH0[locus]/ciHH0[locus]);
    // if(iHS2 > iHS0){
    //     return iHS2;
    // }else{
    //     return 0-iHS0;
    // }
}



/**
 * @brief Calculate the EHH for a single locus as part of IHH routine
 * @param locus The locus 
*/
pair<double, double> IHS::calc_ihh1(int locus){
    int numSnps = hm->mapData->nloci;
    int numHaps = hm->hapData->nhaps;

    if(hm->hapData->unphased){
        double ihh2 = 0;
        double ihh0 = 0;
        double cihh2_upstream = 0;
        double cihh0_upstream  = 0;
        double cihh2_downstream = 0;
        double cihh0_downstream  = 0;
        pair<double, double> ihh2_ihh0_upstream = calc_ehh_unidirection_unphased(locus, false, std::ref(cihh2_upstream), std::ref(cihh0_upstream)); // upstream
        pair<double, double>  ihh2_ihh0_downstream = calc_ehh_unidirection_unphased(locus, true,  std::ref(cihh2_downstream), std::ref(cihh0_downstream)); // downstream
        ihh2 = ihh2_ihh0_upstream.first + ihh2_ihh0_downstream.first;
        ihh0 = ihh2_ihh0_upstream.second + ihh2_ihh0_downstream.second;


        double ihs2 = log10(ihh2/(cihh2_upstream+cihh2_downstream));
        double ihs0 = log10((ihh0/cihh0_upstream+cihh0_downstream));
        
        // if(locus ==0) cout<<"right 0 2 "<<ihh2_ihh0_upstream.second<<" "<<ihh2_ihh0_upstream.first<<endl;
        // if(locus ==0) cout<<"left 0 2 "<<ihh2_ihh0_downstream.second<<" "<<ihh2_ihh0_downstream.first<<endl;
        // if(locus ==0) cout<<"0 not right "<<cihh0_upstream<<" "<<"0 not left" << cihh0_downstream<<endl;
        // if(locus ==0) cout<<"2 not right "<<cihh2_upstream<<" "<<"2 not left" << cihh2_downstream<<endl;
        // if(locus ==0) cout<<"calc "<<ihh2<<" "<<ihh0<<" "<<(cihh2_upstream+cihh2_downstream)<<" " << (cihh0_upstream+cihh0_downstream)<<endl;
        
        return make_pair(ihs2, ihs0);

    }
    

    double ihh1 = 0;
    double ihh0 = 0;

    pair<double, double>  ihh1_ihh0_upstream = calc_ehh_unidirection(locus, false); // upstream
        pair<double, double>  ihh1_ihh0_downstream = calc_ehh_unidirection(locus,  true); // downstream
        ihh1 = ihh1_ihh0_upstream.first + ihh1_ihh0_downstream.first;
        ihh0 = ihh1_ihh0_upstream.second + ihh1_ihh0_downstream.second;
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
    
    }
    return make_pair(ihh1, ihh0);
}

