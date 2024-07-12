#include "xpihh.h"



pthread_mutex_t XPIHH::mutex_log = PTHREAD_MUTEX_INITIALIZER;


void XPIHH::updateEHH_from_split(const unordered_map<unsigned int, vector<unsigned int> > & m, XPIHH_ehh_data* ehhdata){
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


        // if(p.ALT){
        //     del_update = -square_alt(group_count_pooled[old_group_id]) + square_alt(newgroup_size) + square_alt(group_count_pooled[old_group_id] - newgroup_size);
        //     curr_ehh_pooled +=  (del_update  / square_alt(nhaps1+nhaps2)); // if not wagh
        // }else{
        //     if(p.WAGH){

        //     }else{
        //         int del_update = -num_pair(group_count_pooled[old_group_id]) + num_pair(newgroup_size) + num_pair(group_count_pooled[old_group_id] - newgroup_size);
        //         curr_ehh_pooled +=  (del_update / num_pair(nhaps1+nhaps2)); // if not wagh
        //     }   
        // }

            
        
        ehhdata->group_count[old_group_id] -= newgroup_size;
        ehhdata->group_count[ehhdata->totgc] += newgroup_size;  
        
        ehhdata->totgc += 1;
        ehhdata->curr_ehh_before_norm += del_update;
    }
}

/**
 * Calculate EHH in only one direction until cutoff is hit - upstream or downstream
*/
void XPIHH::calc_ehh_unidirection(int locus, unordered_map<unsigned int, vector<unsigned int> > & m, bool downstream){
    XPIHH_ehh_data p1, p2, pooled;
    int numSnps = hm->hapData->nloci; // must be same for both hapData and hapData2
    
    p1.init(hm->hapData->nhaps, hm->hapData->hapEntries[locus].positions);
    p2.init(hm->hapData2->nhaps, hm->hapData2->hapEntries[locus].positions);

    pooled.init(p1.nhaps + p2.nhaps, hm->hapData->hapEntries[locus].positions);
    pooled.init_fix_for_pooled(hm->hapData2->hapEntries[locus].positions, hm->hapData->nhaps);

    p1.initialize_core(p.ALT);
    p2.initialize_core(p.ALT);
    pooled.initialize_core_pooled(p.ALT);

    bool skipLocus = false;

    int i = locus;     //while(true){ // Upstream: for ( int i = locus+1; i<all_positions.size(); i++ )
    while(true){
        if(downstream){
            if (--i < 0) break;
        }else{
            if (++i >= numSnps) break;
        }
        double pooled_ehh = (p.ALT ? pooled.curr_ehh_before_norm / square_alt(pooled.nhaps) : pooled.curr_ehh_before_norm/twice_num_pair(pooled.nhaps));
        //double pooled_ehh = (p.ALT ? p1.curr_ehh_before_norm / square_alt(p1.nhaps) : p1.curr_ehh_before_norm/twice_num_pair(p1.nhaps));
        
        // TODO: should i do for WAGH as well?
        if(pooled_ehh <= p.EHH_CUTOFF){
            break;
        }
        // int nextLocus = i + 1;
        // if(downstream){
        //     nextLocus = i - 1;
        // }
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

        double scale, distance;
        if(p.CALC_XPNSL){
            scale = double(p.SCALE_PARAMETER) / geneticDistance(i, downstream);  
            distance = geneticDistance(i, downstream);
        }else{
            scale = double(p.SCALE_PARAMETER) / physicalDistance(i, downstream);   
            distance = physicalDistance(i, downstream);  
        }
        if(distance > p.SCALE_PARAMETER){
            distance /= p.SCALE_PARAMETER;
        }
        
        if(scale > 1) scale = 1;
        

        if(downstream){
                ihh_p1[locus] += (p1.curr_ehh_before_norm + p1.prev_ehh_before_norm) * 0.5 / twice_num_pair(p1.nhaps);
                ihh_p2[locus] += (p2.curr_ehh_before_norm + p2.prev_ehh_before_norm) * 0.5 / twice_num_pair(p2.nhaps);
        }
        
        p1.prev_ehh_before_norm = p1.curr_ehh_before_norm;   
        p2.prev_ehh_before_norm = p2.curr_ehh_before_norm;    
        pooled.prev_ehh_before_norm = pooled.curr_ehh_before_norm;    

        

        // if(distance> max_gap){
        //     gap_skip = true;
        //     break;
        // }
        
        // if(distance > p.SCALE_PARAMETER){
        //     distance /= p.SCALE_PARAMETER;
        // }
        //distance = 1; // for testing
        
        /*
        if(downstream){
            p1.v = hm->hapData->hapEntries[i+1].xors;
            p2.v = hm->hapData->hapEntries[i+1].xors;
            pooled.v = hm->hapData->hapEntries[i+1].xors;
            pooled.v2 = hm->hapData2->hapEntries[i+1].xors;
        }else{
            p1.v = hm->hapData->hapEntries[i].xors;
            p2.v = hm->hapData->hapEntries[i].xors;
            pooled.v = hm->hapData->hapEntries[i].xors;
            pooled.v2 = hm->hapData2->hapEntries[i].xors;
        }
        */
        p1.v = hm->hapData->hapEntries[i].positions;
        p2.v = hm->hapData2->hapEntries[i].positions;
        pooled.v = p1.v;
        pooled.v2 = p2.v;

        // ensure that in boundary we don't do any calculation
        // if(hm->hapData->hapEntries[i].positions.size() < ones_p1.size() && i!=nhaps1-1 ){ 
        //     ones_p1 = hm->hapData->hapEntries[i].positions;
        //     if(ones_p1.size()==0 or ones_p1.size()==nhaps1){ // integrity check
        //         std::cerr<<"ERROR: Monomorphic site should not exist."<<endl;
        //         throw 0;
        //     }
        // }
        
        // main faster algorithm for ehh

        //for p1
        for (const unsigned int& set_bit_pos : p2.v){
            int old_group_id = p2.group_id[set_bit_pos];
            m[old_group_id].push_back(set_bit_pos);         
        }
        updateEHH_from_split(m, &p2);

        if(p.ALT){
            ihh_p2[locus] += 0.5*distance*(p2.curr_ehh_before_norm + p2.prev_ehh_before_norm)/square_alt(p2.nhaps);
        }else if(p.WAGH){
            ihh_p2[locus] += 0.5*distance*(p2.curr_ehh_before_norm + p2.prev_ehh_before_norm)/(twice_num_pair(p2.n_c[0])+twice_num_pair(p2.n_c[1]));
        }else{
            ihh_p2[locus] += 0.5*distance*(p2.curr_ehh_before_norm + p2.prev_ehh_before_norm)/twice_num_pair(p2.nhaps);
        }
        p2.prev_ehh_before_norm = p2.curr_ehh_before_norm;
        m.clear();

        //for p2
        for (const unsigned int& set_bit_pos : p1.v){
            int old_group_id = p1.group_id[set_bit_pos];
            m[old_group_id].push_back(set_bit_pos);         
        }
        updateEHH_from_split(m, &p1);
        if(p.ALT){
            ihh_p1[locus] += 0.5*distance*(p1.curr_ehh_before_norm + p1.prev_ehh_before_norm)/square_alt(p1.nhaps);
        }else if(p.WAGH){
            ihh_p1[locus] += 0.5*distance*(p1.curr_ehh_before_norm + p1.prev_ehh_before_norm)/(twice_num_pair(p1.n_c[0])+twice_num_pair(p1.n_c[1]));
        }else{
            ihh_p1[locus] += 0.5*distance*(p1.curr_ehh_before_norm + p1.prev_ehh_before_norm)/twice_num_pair(p1.nhaps);
        }
        p1.prev_ehh_before_norm = p1.curr_ehh_before_norm;
        //do not clear m here, as we reuse it in pooled
        
        //for pooled
        for (const unsigned int& set_bit_pos : p2.v){
            int old_group_id = pooled.group_id[set_bit_pos + pooled.nhaps_p1];
            m[old_group_id].push_back(set_bit_pos + pooled.nhaps_p1);         
        }
        updateEHH_from_split(m, &pooled);
        pooled.prev_ehh_before_norm = pooled.curr_ehh_before_norm;
        m.clear(); 
        

        // check if current locus is beyond 1Mb
        if(!p.CALC_XPNSL && physicalDistance_from_core(i,locus, downstream) >= p.MAX_EXTEND) break;
        if(p.CALC_XPNSL && abs(i-locus) >= p.MAX_EXTEND_NSL) break; //g(xi−1, xi) = 1.

        if(skipLocus){
            ihh_p1[locus] = MISSING;
            ihh_p2[locus] = MISSING;
            skipLocus = false;
            continue;
        }
    }
}


void XPIHH::main()
{
    int nloci = hm->mapData->nloci;
    ihh_p1 = new double[nloci];
    ihh_p2 = new double[nloci];

    //barInit(*bar, hm->mapData->nloci, 78); //REWRITE

    if (p.CALC_XPNSL){
        for (int i = 0; i < hm->mapData->nloci; i++){
            hm->mapData->mapEntries[i].geneticPos = i;
        }
    }

    if (p.CALC_XP) std::cerr << "Starting XP-EHH calculations.\n";
    if (p.CALC_XPNSL) std::cerr << "Starting XP-nSL calculations.\n";
        
    std::unordered_map<unsigned int, std::vector<unsigned int> > map_per_thread[numThreads];

    bool openmp_enabled = false; // two different ways to parallelize: first block does pthread, second block does openmp
    if (!openmp_enabled)
    {
        //int total_calc_to_be_done = numSnps;
        std::thread *myThreads = new std::thread[numThreads];
        for (int i = 0; i < numThreads; i++)
        {
            myThreads[i] = std::thread(thread_xpihh, i, std::ref(map_per_thread[i]), this);
        }
        for (int i = 0; i < numThreads; i++)
        {
            myThreads[i].join(); // Join will block our main thread, and so the program won't exit until all finish
        }
        delete[] myThreads;

        (*flog)<<("All threads finished. Now printing xpihh...\n");
    }
    
    std::cerr << "\nFinished.\n";

    if (p.CALC_XP) (*fout) << "id\tpos\tgpos\tp1\tihh1\tp2\tihh2\txpehh\n";
    if (p.CALC_XPNSL) (*fout) << "id\tpos\tgpos\tp1\tsL1\tp2\tsL2\txpnsl\n";
    for (int i = 0; i < hm->mapData->nloci; i++)
    {
        if (ihh_p1[i] != MISSING && ihh_p2[i] != MISSING && ihh_p1[i] != 0 && ihh_p2[i] != 0)
        {
            (*fout) << hm->mapData->mapEntries[i].locusName << "\t"
                    << hm->mapData->mapEntries[i].physicalPos << "\t"
                    << hm->mapData->mapEntries[i].geneticPos << "\t"
                    << hm->hapData->calcFreq(i) << "\t"  //<< freq1[i] << "\t"
                    << ihh_p1[i] << "\t"
                    << hm->hapData2->calcFreq(i) << "\t"  //<< freq2[i] << "\t"
                    << ihh_p2[i] << "\t";
            (*fout) << log10(ihh_p1[i] / ihh_p2[i]) << endl;
        }
    }

    // hm->hapData->releaseHapData();
    // hm->hapData2->releaseHapData();

    delete[] ihh_p1;
    delete[] ihh_p2;
}

/**
 * populate ihh_p1 and ihh_p2 at the end with correct values
*/
void XPIHH::calc_xpihh(int locus, unordered_map<unsigned int, vector<unsigned int> >& m)
{
    ihh_p1[locus] = 0;
    ihh_p2[locus] = 0;
    calc_ehh_unidirection(locus,m, false);
    calc_ehh_unidirection(locus,m,true);   
}

void XPIHH::thread_xpihh(int tid, unordered_map<unsigned int, vector<unsigned int> >& m, XPIHH* obj){
    int numSnps = obj->hm->hapData->nloci;
    int elem_per_block = floor(numSnps/obj->numThreads);
    int start = tid*elem_per_block ;
    int end = start + elem_per_block  ;
    if(tid == obj->numThreads-1 ){
        end = numSnps;
    }

    // int step = (numSnps  / obj->numThreads) / (obj->bar->totalTicks);
    // if (step == 0) step = 1;

    
    //if total 20 tasks: and 4 threads: t0: 0, 4, 8, 12, 16: t1: 1, 5, 9, 13, 17  
    //for (int locus = tid; locus < hm->mapData->nloci; locus += numThreads)
    //#pragma omp parallel 


    for(int locus = start; locus< end; locus++){
        //if (locus % step == 0) advanceBar(*(obj->bar), double(step));
        obj->calc_xpihh(locus, m);
    }

    pthread_mutex_lock(&mutex_log);
    (*(obj->flog))<<("finishing thread # "+to_string(tid)+" at "+to_string(SelscanStats::readTimer())+"\n");
    pthread_mutex_unlock(&mutex_log);
}

