#include "xpihh.h"

pthread_mutex_t XPIHH::mutex_log = PTHREAD_MUTEX_INITIALIZER;
/**
 * Calculate EHH in only one direction until cutoff is hit - upstream or downstream
*/
void XPIHH::calc_ehh_unidirection_xpihh(int locus, unordered_map<unsigned int, vector<unsigned int> > & m, bool downstream){
    HapData& hapData = hm.hapData;
    HapData& hapData2 = hm.hapData2;
    int nhaps1 = hapData.nhaps;
    int nhaps2 = hapData2.nhaps;

    //uint64_t ehh_before_norm = 0;
    //uint64_t curr_ehh_before_norm = 0;

    int group_id_pooled[nhaps1+nhaps2];
    int group_count_pooled[nhaps1+nhaps2];
    int totgc_pooled=0;

    int ancestralCount_p1 = 0;
    int derivedCount_p1 = 0;
    int hetCount_p1 = 0;

    int group_count_p1[nhaps1];
    int group_id_p1[nhaps1];
    bool isDerived_p1[nhaps1];
    bool isAncestral_p1[nhaps1];
    int totgc_p1=0;

    int ancestralCount_p2 = 0;
    int derivedCount_p2 = 0;
    int hetCount_p2 = 0;

    int group_count_p2[nhaps2];
    int group_id_p2[nhaps2];
    bool isDerived_p2[nhaps2];
    bool isAncestral_p2[nhaps2];
    int totgc_p2=0;


    //will be vectorized with compile time flags
    for(int i = 0; i<nhaps1; i++){
        group_count_p1[i] = 0;
        group_id_p1[i] = 0;
        isDerived_p1[i] = false;
        isAncestral_p1[i] = false;
    }

    for(int i = 0; i<nhaps2; i++){
        group_count_p2[i] = 0;
        group_id_p2[i] = 0;
        isDerived_p2[i] = false;
        isAncestral_p2[i] = false;
    }

    for(int i = 0; i<nhaps1+nhaps2; i++){
        group_count_pooled[i] = 0;
        group_id_pooled[i] = 0;
    }

    vector<unsigned int> ones_p1 = hm.hapData.hapEntries[locus].positions;
    vector<unsigned int> ones_p2 = hm.hapData2.hapEntries[locus].positions;

    //TODO
    //vector<unsigned int> twos_p1 = hm.hapData.hapEntries[locus].positions2;
    //vector<unsigned int> twos_p2 = hm.hapData2.hapEntries[locus].positions2;

    
    //unordered_set<unsigned int> v = hm.all_positions[locus];

    //assert : cant be monomorphic
    

    double curr_ehh_p1;
    double prev_ehh_p1;
    double curr_ehh_p2;
    double prev_ehh_p2;
    double curr_ehh_pooled;
    double prev_ehh_pooled;
    double derivedCountPooled, hetCountPooled, ancestralCountPooled;
    bool skipLocus = false;

    

    int i = locus;  
    //while(true){ // Upstream: for ( int i = locus+1; i<all_positions.size(); i++ )
    while(curr_ehh_pooled > p.EHH_CUTOFF){
        int nextLocus = i + 1;
        if(downstream){
            nextLocus = i - 1;
        }
        if (nextLocus < 0 || nextLocus >= hm.mapData.nloci)
        {
            pthread_mutex_lock(&mutex_log);
            (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << p.EHH_CUTOFF
                    << ". ";
            if (!p.TRUNC){
                skipLocus = true;
                (*flog) << "Skipping calculation at position " << hm.mapData.mapEntries[locus].physicalPos << " id: " << hm.mapData.mapEntries[locus].locusName;
            }
            (*flog) << "\n";
            pthread_mutex_unlock(&mutex_log);
            break;
        }
        else if ( physicalDistance(nextLocus, downstream) > p.MAX_GAP )  
        {
            pthread_mutex_lock(&mutex_log);
            (*flog) << "WARNING: Reached a gap of " << physicalDistance(nextLocus, downstream) << "bp > " << p.MAX_GAP 
            << "bp. Skipping calculation at position " <<  hm.mapData.mapEntries[locus].physicalPos << " id: " <<  hm.mapData.mapEntries[locus].locusName << "\n";
            pthread_mutex_unlock(&mutex_log);
            skipLocus = true;
            break;
        }

        double scale = double(p.SCALE_PARAMETER) / physicalDistance(i, downstream);
        if(scale > 1) scale = 1;

        if(i==locus){
            if(hm.hapData.hapEntries[locus].flipped){
                group_count_p1[1] = ones_p1.size();
                group_count_p1[0] = nhaps1 - ones_p1.size();
                ancestralCount_p1 = ones_p1.size();
                derivedCount_p1 = nhaps1 - ones_p1.size();

                for (int set_bit_pos : ones_p1){
                    isAncestral_p1[set_bit_pos] = true;
                    group_id_p1[set_bit_pos] = 1;
                }
                //TODO
            }else{
                group_count_p1[1] = ones_p1.size();
                group_count_p1[0] = nhaps1 - ones_p1.size();

                derivedCount_p1 = ones_p1.size();
                ancestralCount_p1 = nhaps1 - ones_p1.size();

                for (int set_bit_pos : ones_p1){
                    isDerived_p1[set_bit_pos] = true;
                    group_id_p1[set_bit_pos] = 1;
                }
            }        
            
            totgc_p1+=2;
            //ehh0_before_norm = twice_num_pair(n_c0);
            //ehh1_before_norm = twice_num_pair(n_c1);



            // if(downstream){
            //     // if(!calc_all)
            //     //     out_ehh<<"Iter "<<0<<": EHH1["<<locus<<","<<locus<<"]="<<1<<" "<<1<<endl;

            //     if(twice_num_pair(n_c1)!=0){
            //         iHH1[locus] += (curr_ehh1_before_norm + ehh1_before_norm) * 0.5 / twice_num_pair(n_c1);
            //         //cout<<"Summing "<<1.0*curr_ehh1_before_norm/n_c1_squared_minus<<" "<<1.0*ehh1_before_norm/n_c1_squared_minus<<endl;
            //     }
            //     if(twice_num_pair(n_c0)!=0){
            //         iHH0[locus] += (curr_ehh0_before_norm + ehh0_before_norm) * 0.5 / twice_num_pair(n_c0);
            //     }
            // }
            

            //when calculating xp-ehh, ehh does not necessarily start at 1
            if (p.ALT)
            {
                double fD = double(derivedCount_p1) / double(nhaps1);
                double fA = double(ancestralCount_p1) / double(nhaps1);
                double fH = double(hetCount_p1) / double(nhaps1);

                curr_ehh_p1 = fD * fD + fA * fA + fH * fH;
                prev_ehh_p1 = curr_ehh_p1;

                fD = double(derivedCount_p2) / double(nhaps2);
                fA = double(ancestralCount_p2) / double(nhaps2);
                fH = double(hetCount_p2) / double(nhaps2);
                curr_ehh_p2 = fD * fD + fA * fA + fH * fH;
                prev_ehh_p2 = curr_ehh_p2;

                fD = double(derivedCountPooled) / double(nhaps1 + nhaps2);
                fA = double(ancestralCountPooled) / double(nhaps1 + nhaps2);
                fH = double(hetCountPooled) / double(nhaps1 + nhaps2);
                curr_ehh_pooled = fD * fD + fA * fA + fH * fH;
                prev_ehh_pooled = curr_ehh_pooled;
            }
            else
            {
                if (p.WAGH)
                {
                    curr_ehh_p1 = (derivedCount_p1 > 1) ? nCk(derivedCount_p1,2) / (nCk(derivedCount_p1,2)+nCk(nhaps1-derivedCount_p1,2)) : 0;
                    curr_ehh_p1 += (nhaps1 - derivedCount_p1 > 1) ? nCk(nhaps1-derivedCount_p1,2) / (nCk(derivedCount_p1,2)+nCk(nhaps1-derivedCount_p1,2)) : 0;
                    prev_ehh_p1 = curr_ehh_p1;

                    curr_ehh_p2 = (derivedCount_p2 > 1) ? nCk(derivedCount_p2, 2) / (nCk(derivedCount_p2,2)+nCk(nhaps2-derivedCount_p2,2)) : 0;
                    curr_ehh_p2 += (nhaps2 - derivedCount_p2 > 1) ? nCk(nhaps2 - derivedCount_p2, 2) / (nCk(derivedCount_p2,2)+nCk(nhaps2-derivedCount_p2,2)) : 0;
                    prev_ehh_p2 = curr_ehh_p2;

                }
                else
                {
                    curr_ehh_p1 = (derivedCount_p1 > 1) ? nCk(derivedCount_p1, 2) / nCk(nhaps1, 2) : 0;
                    curr_ehh_p1 += (ancestralCount_p1 > 1) ? nCk(ancestralCount_p1, 2) / nCk(nhaps1, 2) : 0;
                    curr_ehh_p1 += (hetCount_p1 > 1) ? nCk(hetCount_p1, 2) / nCk(nhaps1, 2) : 0;

                    prev_ehh_p1 = curr_ehh_p1;

                    curr_ehh_p2 = (derivedCount_p2 > 1) ? nCk(derivedCount_p2, 2) / nCk(nhaps2, 2) : 0;
                    curr_ehh_p2 += (ancestralCount_p2 > 1) ? nCk(ancestralCount_p2, 2) / nCk(nhaps2, 2) : 0;
                    curr_ehh_p2 += (hetCount_p2 > 1) ? nCk(hetCount_p2, 2) / nCk(nhaps2, 2) : 0;
                    prev_ehh_p2 = curr_ehh_p2;        

                }

                curr_ehh_pooled = (derivedCountPooled > 1) ? nCk(derivedCountPooled, 2) / nCk(nhaps1 + nhaps2, 2) : 0;
                curr_ehh_pooled += (ancestralCountPooled > 1) ? nCk(ancestralCountPooled, 2) / nCk(nhaps1 + nhaps2, 2) : 0;
                curr_ehh_pooled += (hetCountPooled > 1) ? nCk(hetCountPooled, 2) / nCk(nhaps1 + nhaps2, 2) : 0;
                prev_ehh_pooled = curr_ehh_pooled;    
            }
        
        }

        if(downstream){
            i--;
        }else{
            i++;
        }
        
        double &distance = scale;

        // if(distance> max_gap){
        //     gap_skip = true;
        //     break;
        // }
        
        // if(distance > p.SCALE_PARAMETER){
        //     distance /= p.SCALE_PARAMETER;
        // }
        //distance = 1; // for testing
        
        if(downstream){
            ones_p1 = hm.hapData.hapEntries[i+1].xors;
            ones_p2 = hm.hapData2.hapEntries[i+1].xors;
        }else{
            ones_p1 = hm.hapData.hapEntries[i].xors;
            ones_p2 = hm.hapData2.hapEntries[i].xors;
        }

        // ensure that in boundary we don't do any calculation
        // if(hm.hapData.hapEntries[i].positions.size() < ones_p1.size() && i!=nhaps1-1 ){ 
        //     ones_p1 = hm.hapData.hapEntries[i].positions;
        //     if(ones_p1.size()==0 or ones_p1.size()==nhaps1){ // integrity check
        //         std::cerr<<"ERROR: Monomorphic site should not exist."<<endl;
        //         throw 0;
        //     }
        // }
        
        // main faster algorithm for ehh
        for (const unsigned int& set_bit_pos : ones_p1){
            int old_group_id = group_id_p1[set_bit_pos];
            m[old_group_id].push_back(set_bit_pos);
        }

        for (const auto &ele : m) {
            int old_group_id = ele.first;
            int newgroup_size = ele.second.size() ;
            if(group_count_p1[old_group_id] == newgroup_size || newgroup_size == 0){
                continue;
            }
            for(int v: ele.second){
                group_id_p1[v] = totgc_p1;
            }
            //int del_update;
            // = -twice_num_pair(group_count_p1[old_group_id]) + twice_num_pair(newgroup_size) + twice_num_pair(group_count_p1[old_group_id] - newgroup_size);
            
            if(p.ALT){
                int del_update = -square_alt(group_count_p1[old_group_id]) + square_alt(newgroup_size) + square_alt(group_count_p1[old_group_id] - newgroup_size);
                curr_ehh_p1 +=  (del_update  / square_alt(nhaps1)); // if not wagh
            }else{
                if(p.WAGH){

                }else{
                    int del_update = -num_pair(group_count_p1[old_group_id]) + num_pair(newgroup_size) + num_pair(group_count_p1[old_group_id] - newgroup_size);
                    curr_ehh_p1 +=  (del_update / num_pair(nhaps1)); // if not wagh
                }   
            }

            group_count_p1[old_group_id] -= newgroup_size;
            group_count_p1[totgc_p1] += newgroup_size;
            totgc_p1+=1;
        }
        ihh_p1[locus] += 0.5*scale*(geneticDistance(locus, downstream))*(curr_ehh_p1 + prev_ehh_p1);
        prev_ehh_p1 = curr_ehh_p1;
        m.clear(); 


        // main faster algorithm for ehh
        for (const unsigned int& set_bit_pos : ones_p2){
            int old_group_id = group_id_p2[set_bit_pos];
            m[old_group_id].push_back(set_bit_pos);
        }

        for (const auto &ele : m) {
            int old_group_id = ele.first;
            int newgroup_size = ele.second.size() ;
            if(group_count_p2[old_group_id] == newgroup_size || newgroup_size == 0){
                continue;
            }
            for(int v: ele.second){
                group_id_p2[v] = totgc_p2;
            }
            //int del_update;
            // = -twice_num_pair(group_count_p1[old_group_id]) + twice_num_pair(newgroup_size) + twice_num_pair(group_count_p1[old_group_id] - newgroup_size);
            
            if(p.ALT){
                int del_update = -square_alt(group_count_p2[old_group_id]) + square_alt(newgroup_size) + square_alt(group_count_p2[old_group_id] - newgroup_size);
                curr_ehh_p2 +=  (del_update  / square_alt(nhaps2)); // if not wagh
            }else{
                if(p.WAGH){

                }else{
                    int del_update = -num_pair(group_count_p2[old_group_id]) + num_pair(newgroup_size) + num_pair(group_count_p2[old_group_id] - newgroup_size);
                    curr_ehh_p2 +=  (del_update / num_pair(nhaps2)); // if not wagh
                }   
            }

            group_count_p1[old_group_id] -= newgroup_size;
            group_count_p1[totgc_p1] += newgroup_size;
            totgc_p1+=1;
        }
        ihh_p2[locus] += 0.5*scale*(geneticDistance(locus, downstream))*(curr_ehh_p2 + prev_ehh_p2);
        prev_ehh_p2 = curr_ehh_p2;
        m.clear(); 

        // main faster algorithm for ehh
        for (const unsigned int& set_bit_pos : ones_p1){
            int old_group_id = group_id_pooled[set_bit_pos];
            m[old_group_id].push_back(set_bit_pos);
        }
        for (const unsigned int& set_bit_pos : ones_p2){
            int old_group_id = group_id_pooled[set_bit_pos];
            m[old_group_id].push_back(set_bit_pos);
        }

        for (const auto &ele : m) {
            int old_group_id = ele.first;
            int newgroup_size = ele.second.size() ;
            if(group_count_pooled[old_group_id] == newgroup_size || newgroup_size == 0){
                continue;
            }
            for(int v: ele.second){
                group_id_pooled[v] = totgc_pooled;
            }
            //int del_update;
            // = -twice_num_pair(group_count_p1[old_group_id]) + twice_num_pair(newgroup_size) + twice_num_pair(group_count_p1[old_group_id] - newgroup_size);
            
            if(p.ALT){
                int del_update = -square_alt(group_count_pooled[old_group_id]) + square_alt(newgroup_size) + square_alt(group_count_pooled[old_group_id] - newgroup_size);
                curr_ehh_pooled +=  (del_update  / square_alt(nhaps1+nhaps2)); // if not wagh
            }else{
                if(p.WAGH){

                }else{
                    int del_update = -num_pair(group_count_pooled[old_group_id]) + num_pair(newgroup_size) + num_pair(group_count_pooled[old_group_id] - newgroup_size);
                    curr_ehh_pooled +=  (del_update / num_pair(nhaps1+nhaps2)); // if not wagh
                }   
            }

            group_count_p1[old_group_id] -= newgroup_size;
            group_count_p1[totgc_p1] += newgroup_size;
            totgc_pooled+=1;
        }
        prev_ehh_pooled = curr_ehh_pooled;
        m.clear(); 

        // check if current locus is beyond 1Mb
        if(!p.CALC_XPNSL && physicalDistance(locus, downstream) >= p.MAX_EXTEND) break;
        if(p.CALC_XPNSL && geneticDistance(locus, downstream) >= p.MAX_EXTEND) break;

        if(skipLocus){
            ihh_p1[locus] = MISSING;
            ihh_p2[locus] = MISSING;
            skipLocus = false;
            continue;
        }
        
        m.clear();
    }
}


void XPIHH::xpihh_main()
{
    int nloci = hm.mapData.nloci;
    ihh_p1 = new double[nloci];
    ihh_p2 = new double[nloci];

    //barInit(*bar, hm.mapData.nloci, 78); //REWRITE

    if (p.CALC_XPNSL){
        for (int i = 0; i < hm.mapData.nloci; i++){
            hm.mapData.mapEntries[i].geneticPos = i;
        }
    }

    if (p.CALC_XP) std::cerr << "Starting XP-EHH calculations.\n";
    if (p.CALC_XPNSL) std::cerr << "Starting XP-nSL calculations.\n";
        

    std::unordered_map<unsigned int, std::vector<unsigned int> > map_per_thread[numThreads];
    std::unordered_map<unsigned int, std::vector<unsigned int> > mapd_per_thread[numThreads];
    bool openmp_enabled = false;
    // two different ways to parallelize: first block does pthread, second block does openmp
    if (!openmp_enabled)
    {
        //int total_calc_to_be_done = numSnps;
        std::thread *myThreads = new std::thread[numThreads];
        for (int i = 0; i < numThreads; i++)
        {
            myThreads[i] = std::thread(thread_xpihh, i, std::ref(map_per_thread[i]),  std::ref(mapd_per_thread[i]), this);
        }
        for (int i = 0; i < numThreads; i++)
        {
            myThreads[i].join(); // Join will block our main thread, and so the program won't exit until all finish
        }
        delete[] myThreads;


        //Logger::write("all threads finished. now calculating ihh...\n");
        (*flog)<<("all threads finished. now printing xpihh...\n");
    }
    
    
    
    if(p.LOW_MEM){
        hm.hapData.releaseHapData_bitset();
        hm.hapData2.releaseHapData_bitset();
    }else{
        hm.hapData.releaseHapData();
        hm.hapData2.releaseHapData();
    }
    
    
    std::cerr << "\nFinished.\n";

    if (p.CALC_XP) (*fout) << "id\tpos\tgpos\tp1\tihh1\tp2\tihh2\txpehh\n";
    if (p.CALC_XPNSL) (*fout) << "id\tpos\tgpos\tp1\tsL1\tp2\tsL2\txpnsl\n";
    for (int i = 0; i < hm.mapData.nloci; i++)
    {
        if (ihh_p1[i] != MISSING && ihh_p2[i] != MISSING && ihh_p1[i] != 0 && ihh_p2[i] != 0)
        {
            (*fout) << hm.mapData.mapEntries[i].locusName << "\t"
                    << hm.mapData.mapEntries[i].physicalPos << "\t"
                    << hm.mapData.mapEntries[i].geneticPos << "\t"
                    << hm.hapData.calcFreq(i) << "\t"  //<< freq1[i] << "\t"
                    << ihh_p1[i] << "\t"
                    << hm.hapData2.calcFreq(i) << "\t"  //<< freq2[i] << "\t"
                    << ihh_p2[i] << "\t";
            (*fout) << log10(ihh_p1[i] / ihh_p2[i]) << endl;
        }
    }

    delete[] ihh_p1;
    delete[] ihh_p2;
}

/**
 * populate ihh_p1 and ihh_p2 at the end with correct values
*/
void XPIHH::calc_xpihh(int locus)
{
    unordered_map<unsigned int, vector<unsigned int> > m;
    ihh_p1[locus] = 0;
    ihh_p2[locus] = 0;
    calc_ehh_unidirection_xpihh(locus,m, false);
    calc_ehh_unidirection_xpihh(locus,m,true);   
}

void XPIHH::thread_xpihh(int tid, unordered_map<unsigned int, vector<unsigned int> >& m, unordered_map<unsigned int, vector<unsigned int> >& md, XPIHH* obj){
    int numSnps = obj->hm.hapData.nloci;
    int elem_per_block = floor(numSnps/obj->numThreads);
    int start = tid*elem_per_block ;
    int end = start + elem_per_block  ;
    if(tid == obj->numThreads-1 ){
        end = numSnps;
    }

    // int step = (numSnps  / obj->numThreads) / (obj->bar->totalTicks);
    // if (step == 0) step = 1;

    
    //if total 20 tasks: and 4 threads: t0: 0, 4, 8, 12, 16: t1: 1, 5, 9, 13, 17  
    //for (int locus = tid; locus < hm.mapData.nloci; locus += numThreads)
    //#pragma omp parallel 
    for(int locus = start; locus< end; locus++){
        //if (locus % step == 0) advanceBar(*(obj->bar), double(step));
        obj->calc_xpihh(locus);
    }

    pthread_mutex_lock(&mutex_log);
    (*(obj->flog))<<("finishing thread # "+to_string(tid)+" at "+to_string(MainTools::readTimer())+"\n");
    pthread_mutex_unlock(&mutex_log);
}

