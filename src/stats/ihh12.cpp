#include "ihh12.h"



pthread_mutex_t IHH12::mutex_log = PTHREAD_MUTEX_INITIALIZER;

// void IHH12::findMaxK(int* arr, int n, int K) {
//     /// Function to find the Kth largest element in the array
//     // using binary search
//     int high = 0; //kthlargest
//     for(int k = 0; k<K; k++){
//     {
//         int low = INT_MAX, high = INT_MIN;

//         // Find the minimum and maximum elements in the array
//         for (int i = 0; i < n; i++) {
//             low = min(low, arr[i]);
//             high = max(high, arr[i]);
//         }

//         // Perform binary search on the range of elements
//         // between low and high
//         while (low <= high) {
//             int mid = low + (high - low) / 2;
//             int count = 0;

//             // Count the number of elements greater than mid in
//             // the array
//             for (int i = 0; i < n; i++) {
//                 if (arr[i] > mid) {
//                     count++;
//                 }
//             }

//             // If there are at least K elements greater than
//             // mid, search the right half
//             if (count >= k) {
//                 low = mid + 1;
//             }
//             // Otherwise, search the left half
//             else {
//                 high = mid - 1;
//             }
//         }

//         // Return the Kth largest element, high will be the
//     }
//     // Print the K largest elements in decreasing order
//     for (int i = 0; i < n; i++) {
//         if (arr[i] >= high) {
//             cout << arr[i] << " ";
//         }
//     }
//     cout << endl;

//      // double firstFreq = (top1 > 1) ? twice_num_pair(top1) : 0;
//     // double secondFreq =(top2 > 1) ?  twice_num_pair(top2): 0;
//     // double comboFreq = ((top1 + top2) > 1) ? twice_num_pair((top1 + top2)) : 0;
//     // double normfac = twice_num_pair(nhaps);
// }

void IHH12::findMaxTwo(int* arr, int n, int &max1, int &max2) {
    // Initialize max1 and max2 to the smallest possible values
    max1 = 0;
    max2 = 0;

    for (int i = 0; i < n; i++) {
        // if (arr[i] > max1) {
        //     max2 = max1;
        //     max1 = arr[i];
        // } else if (arr[i] > max2 && arr[i] != max1) {
        //     max2 = arr[i];
        // }
        if (arr[i]  > max1)
        {
            max2 = max1;
            max1 = arr[i] ;
        }
        else if (arr[i]  > max2)
        {
            max2 = arr[i] ;
        }

    }
}

void IHH12::updateEHH_from_split(const unordered_map<int, vector<int> > & m, IHH12_ehh_data* ehhdata){
    double sum_del_update = 0;
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

        // for(int g = 0; g<ehhdata->totgc ; g++){
        //     if(ehhdata->group_count[g] == 0){
        //         ehhdata->group_count[g] = ehhdata->group_count[ehhdata->totgc-1];
        //         ehhdata->group_count[ehhdata->totgc-1] = 0;
        //         ehhdata->totgc -= 1;
        //         break;
        //     }
        //}
        sum_del_update += del_update;
    }
    ehhdata->curr_ehh_before_norm += sum_del_update;
    int top1, top2;
    findMaxTwo(ehhdata->group_count, ehhdata->totgc, top1, top2);

    double firstFreq = (top1 > 1) ? twice_num_pair(top1) : 0;
    double secondFreq =(top2 > 1) ?  twice_num_pair(top2): 0;
    double comboFreq = ((top1 + top2) > 1) ? twice_num_pair((top1 + top2)) : 0;
    double normfac = twice_num_pair(ehhdata->nhaps);

    ehhdata->curr_ehh12_before_norm = ehhdata->curr_ehh_before_norm  - firstFreq - secondFreq + comboFreq;
    //ehhdata->curr_ehh12_before_norm *= normfac;
    //cout<<"t1 t2 "<<top1<<" "<<top2<<" "<<firstFreq/normfac<<" "<<secondFreq/normfac<<" "<< comboFreq/normfac<<endl;
}

/**
 * Calculate EHH in only one direction until cutoff is hit - upstream or downstream
*/
void IHH12::calc_ehh_unidirection(int locus, unordered_map<int, vector<int> > & m, bool downstream){
    IHH12_ehh_data ehhdata;
    int numSnps = hm->hapData->nloci; // must be same for both hapData and hapData2
    bool &skipLocus = hm->mapData->mapEntries[locus].skipLocus;
    if(skipLocus){ // from downstream
        return;
    }

    ehhdata.init(hm->hapData->nhaps, hm->hapData->hapEntries[locus].hapbitset);
    ehhdata.initialize_core(p.ALT);
    ehhdata.prev_ehh_before_norm = ehhdata.curr_ehh_before_norm;
    ehhdata.prev_ehh12_before_norm = ehhdata.curr_ehh12_before_norm;

    int i = locus;
    while(true){
        // int queryPad = p.QWIN; // only for query locu
        // if(physicalDistance_from_core(i, locus, downstream) > queryPad){
        //     break;
        // }

        if(physicalDistance_from_core(i, locus, downstream) >= max_extend){ //check if currentLocus is beyond 1Mb
            break;
        }

        double breaking_ehh = (p.ALT ? ehhdata.curr_ehh_before_norm / square_alt(ehhdata.nhaps) : ehhdata.curr_ehh_before_norm /twice_num_pair(ehhdata.nhaps));
        if(breaking_ehh <= p.EHH_CUTOFF){
            /*DEBUG :
            if(downstream){
                cout<<"breaking ehh down: "<<locus<<" "<<breaking_ehh<<" "<<ihh12[locus]<<endl;
            }else{
                cout<<"breaking ehh up: "<<locus<<" "<<breaking_ehh<<" "<<ihh12[locus]<<endl;
            }
             */
            break;
        }

        bool breakReachedEdge = false;
        breakReachedEdge = downstream? (i == 0) : (i == numSnps-1);
        if(breakReachedEdge){ //nextLocus < 0 || nextLocus >= numSnps to avoid going to negative
            pthread_mutex_lock(&mutex_log);
            (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << p.EHH_CUTOFF << ". ";
            if (p.TRUNC == false){
                skipLocus = true;
                (*flog) << "Skipping calculation at position " << hm->mapData->mapEntries[i].physicalPos << " id: " << hm->mapData->mapEntries[i].locusName;
            }
            (*flog) << "\n";
            pthread_mutex_unlock(&mutex_log);
            break;
        }

        //at this point we ensured that nextLocus (i-1 or i+1) is within bounds
        i = downstream? i-1 : i+1; //proceed to next locus
        int physicalDistance_old = physicalDistance(i,downstream); //double check if i or nextlocus
        if (physicalDistance_old > p.MAX_GAP)
        {
            pthread_mutex_lock(&mutex_log);
            (*flog) << "WARNING: Reached a gap of " << physicalDistance_old << "bp > " << p.MAX_GAP << "bp. Skipping calculation at position " <<  hm->mapData->mapEntries[i].physicalPos << " id: " <<  hm->mapData->mapEntries[i].locusName << "\n";
            pthread_mutex_unlock(&mutex_log);
            skipLocus = true;
            break;
        }

        double scale, distance;
        if(p.CALC_XPNSL || p.CALC_SOFT){
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
        ehhdata.v = hm->hapData->hapEntries[i].hapbitset;

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
        ACTION_ON_ALL_SET_BITS(ehhdata.v, {
            int old_group_id = ehhdata.group_id[set_bit_pos];
            m[old_group_id].push_back(set_bit_pos); 
        });
        
        updateEHH_from_split(m, &ehhdata);
        m.clear();

        //DEBUG:: if(!downstream) cout<<"before ihh12"<<ihh12[locus]<<endl;
        if(p.ALT){
            ihh12[locus] += 0.5 * (ehhdata.curr_ehh12_before_norm + ehhdata.prev_ehh12_before_norm)/square_alt(ehhdata.nhaps) * scale * distance ;
        }else{
            ihh12[locus] += 0.5  * (ehhdata.curr_ehh12_before_norm/twice_num_pair(ehhdata.nhaps) + ehhdata.prev_ehh12_before_norm/twice_num_pair(ehhdata.nhaps)) * scale * distance;

        }
        double current_ehh1 = ehhdata.curr_ehh_before_norm / twice_num_pair(ehhdata.nhaps);
        double previous_ehh1 = ehhdata.prev_ehh_before_norm / twice_num_pair(ehhdata.nhaps);
        double current_ehh12 = ehhdata.curr_ehh12_before_norm  / twice_num_pair(ehhdata.nhaps);
        double previous_ehh12 = ehhdata.prev_ehh12_before_norm / twice_num_pair(ehhdata.nhaps);
        //if(!downstream) cout<<current_ehh1<<" "<<previous_ehh1<<" "<<current_ehh12<<" "<<previous_ehh12<<" "<<ihh12[locus] <<" "<<scale * distance<< endl;

        ehhdata.prev_ehh_before_norm = ehhdata.curr_ehh_before_norm;
        ehhdata.prev_ehh12_before_norm = ehhdata.curr_ehh12_before_norm;

        if (physicalDistance_from_core(i, locus, downstream) >= max_extend) break;
        
    }
    if(skipLocus==true){
        ihh12[locus] = MISSING;
    }
}


void IHH12::main()
{
    if(p.UNPHASED){
        throw std::invalid_argument("Unphased analysis not yet supported for iHH12 calculations.");
    }

    int nloci = hm->mapData->nloci;
    ihh12 = new double[nloci];

    //load params
    // int SCALE_PARAMETER = p->params->getIntFlag(ARG_GAP_SCALE);
    // int MAX_GAP = p->params->getIntFlag(ARG_MAX_GAP);
    // double EHH_CUTOFF = p->params->getDoubleFlag(ARG_CUTOFF);
    // int numThreads = p->params->getIntFlag(ARG_THREAD);
    // bool TRUNC = p->params->getBoolFlag(ARG_TRUNC);
    // int MAX_EXTEND = ( p->params->getIntFlag(ARG_MAX_EXTEND) <= 0 ) ? physicalPos[nloci - 1] - physicalPos[0] : p->params->getIntFlag(ARG_MAX_EXTEND);


    std::cerr << "Starting iHH12 calculations."<<endl;
    std::unordered_map<int, std::vector<int> >* map_per_thread = new std::unordered_map<int, std::vector<int> > [numThreads];

    bool openmp_enabled = false; // two different ways to parallelize: first block does pthread, second block does openmp
    if (!openmp_enabled)
    {
        //int total_calc_to_be_done = numSnps;
        std::thread *myThreads = new std::thread[numThreads];
        for (int i = 0; i < numThreads; i++)
        {
            myThreads[i] = std::thread(thread_main, i, std::ref(map_per_thread[i]), this);
        }
        for (int i = 0; i < numThreads; i++)
        {
            myThreads[i].join(); // Join will block our main thread, and so the program won't exit until all finish
        }
        delete[] myThreads;

        (*flog)<<("All threads finished. Now printing xpihh...\n");
    }

    std::cerr << "\nFinished.\n";

    delete[] map_per_thread; // free memory

    (*fout) << "id\tpos\tp1\tihh12\n";
    for (int i = 0; i < hm->mapData->nloci; i++)
    {
        if (ihh12[i] != MISSING ) // !hm->mapData->mapEntries[i].skipLocus
        {
            (*fout) << hm->mapData->mapEntries[i].locusName << "\t"
                    << hm->mapData->mapEntries[i].physicalPos << "\t"
                    //<< hm->mapData->mapEntries[i].geneticPos << "\t"
                    << hm->hapData->calcFreq(i) << "\t"  //<< freq1[i] << "\t"
                    << ihh12[i]  << endl;
        }
    }

    // hm->hapData->releaseHapData();
    // hm->hapData2->releaseHapData();
    // hm->mapData->releaseMapData();
    delete[] ihh12;
}

/**
 * populate ihh_p1 and ihh_p2 at the end with correct values
*/
void IHH12::calc_stat_at_locus(int locus, unordered_map<int, vector<int> >& m)
{
    //DEBUG: if(locus==239){
    ihh12[locus] = 0;
    calc_ehh_unidirection(locus, m, true);
    // ihh12[locus] = 0;
    calc_ehh_unidirection(locus, m, false);
}

void IHH12::thread_main(int tid, unordered_map<int, vector<int> >& m, IHH12* obj){
    int numSnps = obj->hm->hapData->nloci;
    int elem_per_block = floor(numSnps/obj->numThreads);
    int start = tid*elem_per_block ;
    int end = start + elem_per_block  ;
    if(tid == obj->numThreads-1 ){
        end = numSnps;
    }

    for(int locus = start; locus< end; locus++){
        obj->calc_stat_at_locus(locus, m);
    }

    pthread_mutex_lock(&mutex_log);
    (*(obj->flog))<<("finishing thread # "+to_string(tid)+" at "+to_string(SelscanStats::readTimer())+"\n");
    pthread_mutex_unlock(&mutex_log);
}

