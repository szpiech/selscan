#include "ihh12.h"
#include <cassert>

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
    
    //double normfac = twice_num_pair(ehhdata->nhaps);

    ehhdata->curr_ehh12_before_norm = ehhdata->curr_ehh_before_norm  - firstFreq - secondFreq + comboFreq;
    //ehhdata->curr_ehh12_before_norm *= normfac;
    //cout<<"t1 t2 "<<top1<<" "<<top2<<" "<<firstFreq/normfac<<" "<<secondFreq/normfac<<" "<< comboFreq/normfac<<endl;
}

/**
 * Calculate EHH in only one direction until cutoff is hit - upstream or downstream
*/
double IHH12::calc_ehh_unidirection(int locus, bool downstream){

    double ihh12 = 0;
    unordered_map<int, vector<int> > m;

    IHH12_ehh_data ehhdata;
    int numSnps = hm->hapData->nloci; // must be same for both hapData and hapData2

    ehhdata.init(hm->hapData->nhaps, hm->hapData->hapEntries[locus].hapbitset);
    ehhdata.initialize_core(p.ALT);
    ehhdata.prev_ehh_before_norm = ehhdata.curr_ehh_before_norm;
    ehhdata.prev_ehh12_before_norm = ehhdata.curr_ehh12_before_norm;

    int i = locus;
    int prev_index = locus; // for multi-MAF
    int skipped_due_to_multimaf = 0; // for multi-MAF

    while(true){

        bool edgeBreak = false;
        edgeBreak = (downstream)? (i-1 < 0) : (i+1 >= numSnps);

        if(edgeBreak) {
            {std::lock_guard<std::mutex> lock(mutex_log);
            (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << p.EHH_CUTOFF
                    << ". position: "<< hm->mapData->mapEntries[locus].physicalPos << " id: " << hm->mapData->mapEntries[locus].locusName << endl;
            if (!p.TRUNC){
                (*flog) << "Skipping calculation.";
            }
            (*flog) << endl;
            }//unlock
            if (!p.TRUNC){
                return skipLocusDouble();
            }
            break;
        }

        i = (downstream)? i-1 : i+1; // update i 

        if(p.MULTI_MAF){
            if(hm->hapData->get_maf(i) < p.MAF){
                skipped_due_to_multimaf++;
                continue;
            }
        }

        // if(physicalDistance(i, locus, downstream) >= max_extend){ //check if currentLocus is beyond 1Mb
        //     break;
        // }




        double distance =  geneticDistance(i, prev_index, downstream);
        double scale = double(p.SCALE_PARAMETER) / double(physicalDistance(i, prev_index, downstream) );
        if (scale > 1) scale = 1;
        distance *= scale;
        if (physicalDistance(i, prev_index, downstream) > p.MAX_GAP)
        {
            {
                std::lock_guard<std::mutex> lock(mutex_log);
                (*flog) << "WARNING: Reached a gap of " << physicalDistance(i, prev_index, downstream)
                        << "bp > " << p.MAX_GAP << "bp. Skipping calculation at position " << hm->mapData->mapEntries[locus].physicalPos << " id: " << hm->mapData->mapEntries[locus].locusName << "\n";
            }
            return skipLocusDouble();
        }
        prev_index = i;  // save last index that was not skipped


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

        // if(distance> max_gap){
        //     gap_skip = true;
        //     break;
        // }

        // if(distance > p.SCALE_PARAMETER){
        //     distance /= p.SCALE_PARAMETER;
        // }
        //distance = 1; // for testing

        ehhdata.v = hm->hapData->hapEntries[i].hapbitset;
        ACTION_ON_ALL_SET_BITS(ehhdata.v, {
            int old_group_id = ehhdata.group_id[set_bit_pos];
            m[old_group_id].push_back(set_bit_pos); 
        });
        
        updateEHH_from_split(m, &ehhdata);
        m.clear();

        //DEBUG:: if(!downstream) cout<<"before ihh12"<<ihh12[locus]<<endl;
        if(p.ALT){
            ihh12 += 0.5 * (ehhdata.curr_ehh12_before_norm + ehhdata.prev_ehh12_before_norm)/square_alt(ehhdata.nhaps) * scale * distance ;
        }else{
            ihh12 += 0.5  * (ehhdata.curr_ehh12_before_norm/twice_num_pair(ehhdata.nhaps) + ehhdata.prev_ehh12_before_norm/twice_num_pair(ehhdata.nhaps)) * scale * distance;
        }

        // double current_ehh1 = ehhdata.curr_ehh_before_norm / twice_num_pair(ehhdata.nhaps);
        // double previous_ehh1 = ehhdata.prev_ehh_before_norm / twice_num_pair(ehhdata.nhaps);

        // double current_ehh12 = ehhdata.curr_ehh12_before_norm  / twice_num_pair(ehhdata.nhaps);
        // double previous_ehh12 = ehhdata.prev_ehh12_before_norm / twice_num_pair(ehhdata.nhaps);
        
        //if(!downstream) cout<<current_ehh1<<" "<<previous_ehh1<<" "<<current_ehh12<<" "<<previous_ehh12<<" "<<ihh12[locus] <<" "<<scale * distance<< endl;

        ehhdata.prev_ehh_before_norm = ehhdata.curr_ehh_before_norm;
        ehhdata.prev_ehh12_before_norm = ehhdata.curr_ehh12_before_norm;

        if(ehhdata.totgc == hm->hapData->nhaps) break;
        if (physicalDistance(i, locus, downstream) >= max_extend) break;
        
    }
    return ihh12;
}


void IHH12::main()
{
    init_global_fout("ihh12");

    if(p.UNPHASED){
        HANDLE_ERROR("Unphased analysis not supported for iHH12 calculations.");
    }


    std::cerr << "Starting iHH12 calculations."<<endl;


    ThreadPool pool(p.numThreads);
    std::vector< std::future<double> > results;
    for(int i = 0; i <  hm->mapData->nloci; ++i) {
        results.emplace_back(
            pool.enqueue([i,this] {
                return this->calc_ihh12_at_locus(i);
            })
        );
    }

    int i = 0;
    (*fout) << "id\tpos\tp1\tihh12\n";
    for(auto && result: results){ // this is a blocking call
        double ihh12 = result.get(); 

        if(ihh12 != skipLocusDouble()){
            (*fout) << hm->mapData->mapEntries[i].locusName << "\t"
                    << hm->mapData->mapEntries[i].physicalPos << "\t"
                    //<< hm->mapData->mapEntries[i].geneticPos << "\t"
                    << hm->hapData->calcFreq(i) << "\t"  //<< freq1[i] << "\t"
                    << ihh12  << endl;     
        }

        i++;
        assert(i <= hm->mapData->nloci);
    }
    //close fouts
    fout->close();
    // hm->hapData->releaseHapData();
    // hm->hapData2->releaseHapData();
    // hm->mapData->releaseMapData();

}

/**
 * populate ihh_p1 and ihh_p2 at the end with correct values
*/
double IHH12::calc_ihh12_at_locus(int locus)
{
    if(hm->hapData->get_maf(locus) < p.MAF) { // if core locus has MAF < p.MAF, skip it, useful in --keep-low-freq
    {
            string reason = "MAF < " + std::to_string(p.MAF);
            std::lock_guard<std::mutex> lock(mutex_log);
            (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
            << " (pos " << hm->mapData->mapEntries[locus].physicalPos << ") has " 
            << reason << ". Skipping this locus." << std::endl;
        }
    return skipLocusDouble();
    }


    double ihh12_down = calc_ehh_unidirection(locus, true);
    if(ihh12_down == skipLocusDouble()){
        return skipLocusDouble();
    }

    double ihh12_up = calc_ehh_unidirection(locus, false);
    if(ihh12_up == skipLocusDouble()){
        return skipLocusDouble();
    }

    return ihh12_down + ihh12_up;
}

