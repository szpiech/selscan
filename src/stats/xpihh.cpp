#include "xpihh.h"
#include "../thread_pool.h"

/* Implementation log */
/// LOW_MEM AND HAP
/// NOT LOW MEM AND HAP AND VCF
/// NO UNPHASED

pthread_mutex_t XPIHH::mutex_log = PTHREAD_MUTEX_INITIALIZER;


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

XPIHH_ehh_data::~XPIHH_ehh_data(){
    delete[] group_count;
    delete[] group_id;
}

void XPIHH::updateEHH_from_split(const unordered_map<unsigned int, vector<unsigned int> > & m, XPIHH_ehh_data* ehhdata){
    for (const auto &ele : m) {
        int old_group_id = ele.first;
        int newgroup_size = ele.second.size() ;

        if(ehhdata->group_count[old_group_id] == newgroup_size || newgroup_size == 0){ //if a group becomes empty, we don't increment num_groups, we just reuse that group
            continue;
        }

        for(int v: ele.second){
            ehhdata->group_id[v] = ehhdata->totgc;
            if(v > ehhdata->nhaps || v < 0){
                throw std::runtime_error("set_bit_pos out of bounds");
            }
        }
        
        double del_update = -twice_num_pair(ehhdata->group_count[old_group_id]) + twice_num_pair(newgroup_size) + twice_num_pair(ehhdata->group_count[old_group_id] - newgroup_size);
        if(p.ALT){
            del_update = -square_alt(ehhdata->group_count[old_group_id]) +   square_alt(newgroup_size) + square_alt(ehhdata->group_count[old_group_id] - newgroup_size);
        }
        
        ehhdata->group_count[old_group_id] -= newgroup_size;
        if(old_group_id > ehhdata->nhaps || old_group_id < 0){
            cout<<old_group_id<<" "<<ehhdata->nhaps<<endl;
                throw std::runtime_error("set_bit_pos out of bounds");
            }
        ehhdata->group_count[ehhdata->totgc] += newgroup_size;  
        if(ehhdata->totgc > ehhdata->nhaps || ehhdata->totgc < 0){
                cout<<ehhdata->totgc<<" "<<ehhdata->nhaps<<endl;

                throw std::runtime_error("set_bit_pos out of bounds");
            }
        
        ehhdata->totgc += 1;
        ehhdata->curr_ehh_before_norm += del_update;
        if(ehhdata->totgc > ehhdata->nhaps){
            cerr<<ehhdata->totgc<<" "<<ehhdata->nhaps<<endl;
            throw std::runtime_error("totgc out of bounds");
            exit(1);
        }
    }
}

// ----------------------------------- //

/**
 * Calculate EHH in only one direction until cutoff is hit - upstream or downstream
*/
pair<double, double> XPIHH::calc_ehh_unidirection(int locus,  bool downstream){
    std::unique_ptr<std::unordered_map<unsigned int, std::vector<unsigned int>>> mp(new std::unordered_map< unsigned int, std::vector<unsigned int>>());
    unordered_map<unsigned int, vector<unsigned int> >& m = (* mp);

    int numSnps = hm->hapData->nloci; // must be same for both hapData and hapData2
    if(hm->hapData->nloci != hm->hapData2->nloci){
        throw std::runtime_error("Number of SNPs in both hapData and hapData2 should be same.");
    }

    double ihh_p1 = 0;
    double ihh_p2 = 0;

    XPIHH_ehh_data* p1, *p2, *pooled;
    p1 = new XPIHH_ehh_data();
    p2 = new XPIHH_ehh_data();
    pooled = new XPIHH_ehh_data();
    
    p1->nhaps = hm->hapData->nhaps;
    p2->nhaps = hm->hapData2->nhaps;
    pooled->nhaps = p1->nhaps + p2->nhaps;

    p1->group_count = new int[p1->nhaps];
    p2->group_count = new int[p2->nhaps];
    pooled->group_count = new int[pooled->nhaps];

    p1->group_id = new int[p1->nhaps];
    p2->group_id = new int[p2->nhaps];
    pooled->group_id = new int[pooled->nhaps];

    for(int i = 0; i< p1->nhaps; i++){
        p1->group_count[i] = 0;
        p1->group_id[i] = 0;
        pooled->group_count[i] = 0;
        pooled->group_id[i] = 0;
    }

    for(int i = 0; i< p2->nhaps; i++){
        p2->group_count[i] = 0;
        p2->group_id[i] = 0;
        pooled->group_count[p1->nhaps + i] = 0;
        pooled->group_id[p1->nhaps + i] = 0;
    }

    vector<unsigned int>* v1 = &hm->hapData->hapEntries[locus].positions;
    vector<unsigned int>* v2 = &hm->hapData2->hapEntries[locus].positions;
    MyBitset* vb1 = hm->hapData->hapEntries[locus].hapbitset;
    MyBitset* vb2 = hm->hapData2->hapEntries[locus].hapbitset;

    if(p.LOW_MEM){
        p1->n_c[1] = vb1->num_1s;
        p2->n_c[1] = vb2->num_1s;
    }else{
        p1->n_c[1] = (*v1).size();
        p2->n_c[1] = (*v2).size();
    }
    
    p1->n_c[0] = p1->nhaps - p1->n_c[1];
    p2->n_c[0] = p2->nhaps - p2->n_c[1];

    pooled->n_c[1] = p1->n_c[1] + p2->n_c[1];
    pooled->n_c[0] = pooled->nhaps - pooled->n_c[1];

    if(p.ALT){
        p1->normalizer = square_alt(p1->nhaps);
        p2->normalizer = square_alt(p2->nhaps);
        pooled->normalizer = square_alt(pooled->nhaps);
    }else if(p.WAGH){
        p1->normalizer = twice_num_pair(p1->n_c[0])+twice_num_pair(p1->n_c[1]);
        p2->normalizer = twice_num_pair(p2->n_c[0])+twice_num_pair(p2->n_c[1]);
        pooled->normalizer = twice_num_pair(pooled->n_c[0])+twice_num_pair(pooled->n_c[1]);
    }else{
        p1->normalizer = twice_num_pair(p1->nhaps);
        p2->normalizer = twice_num_pair(p2->nhaps);
        pooled->normalizer = twice_num_pair(pooled->nhaps);
    }

    //init core
    if(p1->n_c[1]==0 || p1->n_c[1]==p1->nhaps){    //monomorphic
        p1->group_count[0] = p1->nhaps;
        p1->totgc+=1;
    }else{
        p1->group_count[1] = p1->n_c[1];
        p1->group_count[0] = p1->n_c[0];
        p1->totgc+=2;

        if(p.LOW_MEM){
            ACTION_ON_ALL_SET_BITS(vb1, {
                p1->group_id[set_bit_pos] = 1;
            });
        }else{
            for (int set_bit_pos : *(v1)){
                p1->group_id[set_bit_pos] = 1;
            }
        }
    }

    if(p2->n_c[1]==0 || p2->n_c[1]==p2->nhaps){    //monomorphic
        p2->group_count[0] = p2->nhaps;
        p2->totgc+=1;
    }else{
        p2->group_count[1] = p2->n_c[1];
        p2->group_count[0] = p2->n_c[0];
        p2->totgc+=2;

        if(p.LOW_MEM){
            ACTION_ON_ALL_SET_BITS(vb2, {
                p2->group_id[set_bit_pos] = 1;
            });
        }else{
            for (int set_bit_pos : (*v2)){
                p2->group_id[set_bit_pos] = 1;
            }
        }
        
    }

    if(pooled->n_c[1]==0 || pooled->n_c[1]==pooled->nhaps){    //monomorphic
        pooled->group_count[0] = pooled->nhaps;
        pooled->totgc+=1;
    }else{
        pooled->group_count[1] = pooled->n_c[1];
        pooled->group_count[0] = pooled->n_c[0];
        pooled->totgc+=2;

        if(p.LOW_MEM){
            ACTION_ON_ALL_SET_BITS(vb1, {
                pooled->group_id[set_bit_pos] = 1;
            });
            ACTION_ON_ALL_SET_BITS(vb2, {
                pooled->group_id[set_bit_pos + p1->nhaps] = 1;
            });
        }else{
            for (int set_bit_pos : *(v1)){
                pooled->group_id[set_bit_pos] = 1;
            }
            for (int set_bit_pos : (*v1)){
                pooled->group_id[set_bit_pos + p1->nhaps] = 1;
            }
        }

    }

        

    if(p.ALT){
        p1->curr_ehh_before_norm = square_alt(p1->n_c[1])+square_alt(p1->n_c[0]);
        p2->curr_ehh_before_norm = square_alt(p2->n_c[1])+square_alt(p2->n_c[0]);
        pooled->curr_ehh_before_norm = square_alt(pooled->n_c[1])+square_alt(pooled->n_c[0]);
    }else{
        p1->curr_ehh_before_norm = twice_num_pair(p1->n_c[1])+twice_num_pair(p1->n_c[0]);
        p2->curr_ehh_before_norm = twice_num_pair(p2->n_c[1])+twice_num_pair(p2->n_c[0]);
        pooled->curr_ehh_before_norm = twice_num_pair(pooled->n_c[1])+twice_num_pair(pooled->n_c[0]);
    }

    p1->prev_ehh_before_norm = p1->curr_ehh_before_norm;
    p2->prev_ehh_before_norm = p2->curr_ehh_before_norm;
    pooled->prev_ehh_before_norm = pooled->curr_ehh_before_norm;
    
    int i = locus;     
    while(true){
        double pooled_ehh = (p.ALT ? pooled->curr_ehh_before_norm / square_alt(pooled->nhaps) : pooled->curr_ehh_before_norm/twice_num_pair(pooled->nhaps));
        
        // TODO: should i do for WAGH as well?
        if(pooled_ehh <= p.EHH_CUTOFF){
            break;
        }

        if ((downstream && i-1<0) || (!downstream && i+1>=numSnps))
        {
            pthread_mutex_lock(&mutex_log);
            (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << p.EHH_CUTOFF
                    << ". ";
            if (!p.TRUNC){
                 hm->mapData->mapEntries[locus].skipLocus = true;
                (*flog) << "Skipping calculation at position " << hm->mapData->mapEntries[locus].physicalPos << " id: " << hm->mapData->mapEntries[locus].locusName;
            }
            (*flog) << "\n";
            pthread_mutex_unlock(&mutex_log);
            break;
        }

        if(downstream){
            if (--i < 0) break;
        }else{
            if (++i >= numSnps) break;
        }

        if (physicalDistance(i,downstream) > p.MAX_GAP)
        {
            pthread_mutex_lock(&mutex_log);
            (*flog) << "WARNING: Reached a gap of " << physicalDistance(i,downstream)
                    << "bp > " << p.MAX_GAP << "bp. Skipping calculation at position " << hm->mapData->mapEntries[locus].physicalPos << " id: " << hm->mapData->mapEntries[locus].locusName << "\n";
            pthread_mutex_unlock(&mutex_log);
            hm->mapData->mapEntries[locus].skipLocus = true;
            break;
        }

        double scale, distance;
        if(p.CALC_XPNSL){
            distance = 1;
        }else{
            distance = physicalDistance(i, downstream);  
        }
        scale = double(p.SCALE_PARAMETER) / physicalDistance(i, downstream);   
        if(scale > 1) scale = 1;
        distance *= scale;
        
        int i_xor = downstream ? i+1 : i;
        size_t xor_size; 
        size_t pos_size;
        vector<unsigned int>* v; 
        MyBitset* vb; 

        if(p.LOW_MEM){
            xor_size = hm->hapData->hapEntries[i_xor].xorbitset->num_1s;
            pos_size = hm->hapData->hapEntries[i].hapbitset->num_1s;
            vb =  (xor_size>pos_size) ? hm->hapData->hapEntries[i].hapbitset : hm->hapData->hapEntries[i_xor].xorbitset;
        }else{
            xor_size = hm->hapData->hapEntries[i_xor].xors.size();
            pos_size = hm->hapData->hapEntries[i].positions.size();
            v =  (xor_size>pos_size) ? &(hm->hapData->hapEntries[i].positions) : &(hm->hapData->hapEntries[i_xor].xors);
        }

        if(p1->totgc != p1->nhaps){
            if(p.LOW_MEM){
                ACTION_ON_ALL_SET_BITS(vb, {
                    int old_group_id = p1->group_id[set_bit_pos];
                    m[old_group_id].push_back(set_bit_pos);      
                });
            }else{
                for (const unsigned int& set_bit_pos : *v){
                    int old_group_id = p1->group_id[set_bit_pos];
                    m[old_group_id].push_back(set_bit_pos);        
                }
            }
            updateEHH_from_split(m, p1);
            ihh_p1 += 0.5*distance*(p1->curr_ehh_before_norm + p1->prev_ehh_before_norm)/p1->normalizer;
    
            m.clear();

            p1->prev_ehh_before_norm = p1->curr_ehh_before_norm;   
        }

        if(p.LOW_MEM){
            xor_size = hm->hapData2->hapEntries[i_xor].xorbitset->num_1s;
            pos_size = hm->hapData2->hapEntries[i].hapbitset->num_1s;
            vb =  (xor_size>pos_size) ? hm->hapData2->hapEntries[i].hapbitset : hm->hapData2->hapEntries[i_xor].xorbitset;
        }else{
            xor_size = hm->hapData2->hapEntries[i_xor].xors.size();
            pos_size = hm->hapData2->hapEntries[i].positions.size();
            v =  (xor_size>pos_size) ? &(hm->hapData2->hapEntries[i].positions) : &(hm->hapData2->hapEntries[i_xor].xors);
        }
        
        if(p2->totgc != p2->nhaps){
            if(p.LOW_MEM){
                ACTION_ON_ALL_SET_BITS(vb, {
                    int old_group_id = p2->group_id[set_bit_pos];
                    m[old_group_id].push_back(set_bit_pos);  
                });
            }else{
                for (const unsigned int& set_bit_pos : *v){
                    int old_group_id = p2->group_id[set_bit_pos];
                    m[old_group_id].push_back(set_bit_pos);         
                }
            }
            
            updateEHH_from_split(m, p2);
            ihh_p2 += 0.5*distance*(p2->curr_ehh_before_norm + p2->prev_ehh_before_norm)/p2->normalizer;
            p2->prev_ehh_before_norm = p2->curr_ehh_before_norm;
            m.clear();
        }


        if(p.LOW_MEM){
            xor_size = hm->hapData2->hapEntries[i_xor].xorbitset->num_1s + hm->hapData->hapEntries[i_xor].xorbitset->num_1s;
            pos_size = hm->hapData2->hapEntries[i].hapbitset->num_1s + hm->hapData->hapEntries[i].hapbitset->num_1s;
        }else{
            xor_size = hm->hapData->hapEntries[i_xor].xors.size() + hm->hapData2->hapEntries[i_xor].xors.size();
            pos_size =  hm->hapData->hapEntries[i].positions.size() + hm->hapData2->hapEntries[i].positions.size();
        }

        if(pooled->totgc != pooled->nhaps){
            if(p.LOW_MEM){
                //NOW POOLED 1
                vb =  (xor_size>pos_size) ? hm->hapData->hapEntries[i].hapbitset : hm->hapData->hapEntries[i_xor].xorbitset;
                ACTION_ON_ALL_SET_BITS(vb, {
                    int old_group_id = pooled->group_id[set_bit_pos];
                    m[old_group_id].push_back(set_bit_pos); 
                });
                //for pooled2
                vb =  (xor_size>pos_size) ? hm->hapData2->hapEntries[i].hapbitset : hm->hapData2->hapEntries[i_xor].xorbitset;
                ACTION_ON_ALL_SET_BITS(vb, {
                    int old_group_id = pooled->group_id[set_bit_pos + p1->nhaps];
                    m[old_group_id].push_back(set_bit_pos + p1->nhaps);  
                });

            }else{
                //NOW POOLED 1
                v =  (xor_size>pos_size) ? &(hm->hapData->hapEntries[i].positions) : &(hm->hapData->hapEntries[i_xor].xors);
                for (const unsigned int& set_bit_pos : *v){
                    int old_group_id = pooled->group_id[set_bit_pos];
                    m[old_group_id].push_back(set_bit_pos); 
                }
                //for pooled2
                v =  (xor_size>pos_size) ? &(hm->hapData2->hapEntries[i].positions) : &(hm->hapData2->hapEntries[i_xor].xors);
                for (const unsigned int& set_bit_pos : *v){
                    int old_group_id = pooled->group_id[set_bit_pos + p1->nhaps];
                    m[old_group_id].push_back(set_bit_pos + p1->nhaps);  
                }
            }
            
            updateEHH_from_split(m, pooled);
            pooled->prev_ehh_before_norm = pooled->curr_ehh_before_norm;
            // no need to update pooled ihh as we are not interested in it
            m.clear(); 
        }

        // check if current locus is beyond 1Mb
        if(!p.CALC_XPNSL && physicalDistance_from_core(i,locus, downstream) >= max_extend) break;
        if(p.CALC_XPNSL && abs(i-locus) >= max_extend) break; //g(xi−1, xi) = 1.
        // if(downstream){
        //     cout<<locus<<":::l "<<i << " "<<pooled->totgc<<" "<<p1->totgc<<" "<<p2->totgc<<"  p"<<pooled->curr_ehh_before_norm/pooled->normalizer<<" "<<p1->curr_ehh_before_norm/p1->normalizer<<" "<<p2->curr_ehh_before_norm/p2->normalizer<<" "<<ihh_p1<<" "<<ihh_p2<<endl;

        // }else{
        //     cout<<locus<<":::r "<<i << " "<<pooled->totgc<<" "<<p1->totgc<<" "<<p2->totgc<<"  p"<<pooled->curr_ehh_before_norm/pooled->normalizer<<" "<<p1->curr_ehh_before_norm/p1->normalizer<<" "<<p2->curr_ehh_before_norm/p2->normalizer<<" "<<ihh_p1<<" "<<ihh_p2<<endl;
        // }
    }

    delete p1;
    delete p2;
    delete pooled;

    return make_pair(ihh_p1, ihh_p2);
}

void XPIHH::main()
{
    int nloci = hm->mapData->nloci;

    if (p.CALC_XPNSL){
        for (int i = 0; i < hm->mapData->nloci; i++){
            hm->mapData->mapEntries[i].geneticPos = i;
        }
    }

    if (p.CALC_XP) std::cerr << "Starting XP-EHH calculations.\n";
    if (p.CALC_XPNSL) std::cerr << "Starting XP-nSL calculations.\n";

    if (p.CALC_XP) (*fout) << "id\tpos\tgpos\tp1\tihh1\tp2\tihh2\txpehh\n";
    if (p.CALC_XPNSL) (*fout) << "id\tpos\tgpos\tp1\tsL1\tp2\tsL2\txpnsl\n";

    if(numThreads==1){
        for (int i = 0; i < nloci; i++)
        {
            pair<double, double> ihh_p1_p2 = calc_xpihh(i);
            double ihh_p1 = ihh_p1_p2.first;
            double ihh_p2 = ihh_p1_p2.second;
            
            if ( hm->mapData->mapEntries[i].skipLocus == false && ihh_p1 != 0 && ihh_p2 != 0)
            {
                (*fout) << hm->mapData->mapEntries[i].locusName << "\t"
                        << hm->mapData->mapEntries[i].physicalPos << "\t"
                        << hm->mapData->mapEntries[i].geneticPos << "\t"
                        << hm->hapData->calcFreq(i) << "\t"  //<< freq1[i] << "\t"
                        << ihh_p1 << "\t"
                        << hm->hapData2->calcFreq(i) << "\t"  //<< freq2[i] << "\t"
                        << ihh_p2 << "\t";
                (*fout) << log10(ihh_p1 / ihh_p2) << endl;
            }
            
        }
    }else{
        ThreadPool pool(p.numThreads);
        std::vector< std::future<pair<double, double> > > results;
        for(int i = 0; i <  hm->mapData->nloci; ++i) {
            results.emplace_back(
                pool.enqueue([i,this] {
                    return this->calc_xpihh(i);
                })
            );
        }
        int locus = 0;
        for(auto && result: results){ // this is a blocking call
            pair<double, double> ihh_p1_p2 = result.get(); 
            double ihh_p1 = ihh_p1_p2.first;
            double ihh_p2 = ihh_p1_p2.second;

            if ( hm->mapData->mapEntries->skipLocus == false && ihh_p1 != 0 && ihh_p2 != 0)
            {
                (*fout) << hm->mapData->mapEntries[locus].locusName << "\t"
                        << hm->mapData->mapEntries[locus].physicalPos << "\t"
                        << hm->mapData->mapEntries[locus].geneticPos << "\t"
                        << hm->hapData->calcFreq(locus) << "\t"  //<< freq1[i] << "\t"
                        << ihh_p1 << "\t"
                        << hm->hapData2->calcFreq(locus) << "\t"  //<< freq2[i] << "\t"
                        << ihh_p2 << "\t";
                (*fout) << log10(ihh_p1 / ihh_p2) << endl;
            }

            locus++;
            if(locus>hm->mapData->nloci){
                throw std::runtime_error("locus out of bounds");
            }
        }
    }
    std::cerr << "\nFinished.\n";
}

/**
 * populate ihh_p1 and ihh_p2 at the end with correct values
*/
pair<double, double> XPIHH::calc_xpihh(int locus)
{
    pair<double, double> right = calc_ehh_unidirection(locus, false);
    pair<double, double> left = calc_ehh_unidirection(locus, true);   
    return make_pair(left.first+right.first, left.second+right.second);
}
