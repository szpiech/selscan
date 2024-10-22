// #include "ehh12.h"


// /**
//  * Calculate EHH in only one direction until cutoff is hit - upstream or downstream
// */
// void EHH12::calc_ehh_unidirection(int locus, unordered_map<int, vector<int> > & m, bool downstream){
//     int numHaps = hm.hapData.nhaps;
//     int numSnps = hm.mapData.nloci;
//     bool unphased = p.UNPHASED;

//     int total_iteration_of_m = 0;
//     double ehh_before_norm = 0;
//     double ehh0_before_norm = 0;
//     double ehh1_before_norm = 0;

//     bool gap_skip = false;

//     uint64_t n_c0=0;
//     uint64_t n_c1=0;

//     uint64_t &ancestralCount = n_c0;
//     uint64_t &derivedCount = n_c1;

//     uint64_t core_n_c0=0;
//     uint64_t core_n_c1=0;

//     int group_count[numHaps];
//     int group_id[numHaps];
//     bool isDerived[numHaps];
//     bool isAncestral[numHaps];

//     //will be vectorized with compile time flags
//     for(int i = 0; i<numHaps; i++){
//         group_count[i] = 0;
//         group_id[i] = 0;
//         isDerived[i] = false;
//         isAncestral[i] = false;
//     }

//     int totgc=0;
//     vector<int> v = hm.hapData.hapEntries[locus].positions;
    
//     if(v.size() == 0 or  v.size() == numHaps){
//         std::cerr<<"ERROR: Monomorphic site should not exist";
//         throw 0;
//     }

//     hm.hapData.hapEntries[locus].flipped = false;
//     if(hm.hapData.hapEntries[locus].flipped){
//         group_count[1] = v.size();
//         group_count[0] = numHaps - v.size();
//         n_c0 = v.size();
//         n_c1 = numHaps - v.size();

//         for (int set_bit_pos : v){
//             isAncestral[set_bit_pos] = true;
//             group_id[set_bit_pos] = 1;
//         }
//     }else{
//         group_count[1] = v.size();
//         group_count[0] = numHaps - v.size();
//         n_c1 = v.size();
//         n_c0 = numHaps - v.size();

//         for (int set_bit_pos : v){
//             isDerived[set_bit_pos] = true;
//             group_id[set_bit_pos] = 1;
//         }
//     }        
//     totgc+=2;
//     ehh0_before_norm = twice_num_pair(n_c0);
//     ehh1_before_norm = twice_num_pair(n_c1);
//     ehh_before_norm = (twice_num_pair(n_c0)+twice_num_pair(n_c1));
//     //twice_num_pair(n_c1+n_c0);

    
//     int i = locus;  
//     while(true){ // Upstream: for ( int i = locus+1; i<all_positions.size(); i++ )
//         if(downstream){
//             if (--i < 0) break;
//             //if (hm.mentries[locus].phyPos - hm.mentries[i].phyPos > max_extend) break;
//         }else{
//             if (++i >= numSnps) break;
//             //if (hm.mentries[i].phyPos -hm.mentries[locus].phyPos > max_extend) break;
//         }
        
        
//         //if(curr_ehh1_before_norm*1.0/n_c1_squared_minus < cutoff and curr_ehh0_before_norm*1.0/n_c0_squared_minus < cutoff){
        
//         ///*
//         if(ehh1_before_norm*1.0/twice_num_pair(n_c1) <= p.EHH_CUTOFF or ehh0_before_norm*1.0/twice_num_pair(n_c0)  <= p.EHH_CUTOFF){   // or cutoff, change for benchmarking against hapbin
//             //std::cout<<"breaking"<<endl;
//             break;
//         }
//         //*/
//         double distance;
        
//         if(downstream){
//             distance = hm.mapData.mapEntries[i+1].physicalPos - hm.mapData.mapEntries[i].physicalPos;
//         }else{
//             distance = hm.mapData.mapEntries[i].physicalPos - hm.mapData.mapEntries[i-1].physicalPos;
//         }

//         // this should not happen as we already did integrity check previously
//         if (distance < 0)
//         {
//             std::cerr << "ERROR: physical position not in ascending order.\n"; 
//             throw 0;
//         }
        
        
//         if(downstream){
//             v = hm.hapData.hapEntries[i+1].xors;
//         }else{
//             v = hm.hapData.hapEntries[i].xors;
//         }

//         if(hm.hapData.hapEntries[i].positions.size() < v.size() && i!=numHaps-1 ){ //  dont do in boundary
//             v = hm.hapData.hapEntries[i].positions;
//         }
        
//         for (const int& set_bit_pos : v){
//             int old_group_id = group_id[set_bit_pos];
//             m[old_group_id].push_back(set_bit_pos);
//         }

//         for (const auto &ele : m) {
            
//             int old_group_id = ele.first;
//             int newgroup_size = ele.second.size() ;
                            
//             total_iteration_of_m += newgroup_size;
                            
//             if(group_count[old_group_id] == newgroup_size || newgroup_size == 0){
//                 continue;
//             }

//             for(int v: ele.second){
//                 group_id[v] = totgc;
//             }
            
//             double del_update = -twice_num_pair(group_count[old_group_id]) + twice_num_pair(newgroup_size) + twice_num_pair(group_count[old_group_id] - newgroup_size);
//             if(p.ALT){
//                 del_update = -square_alt(group_count[old_group_id]) +   square_alt(newgroup_size) + square_alt(group_count[old_group_id] - newgroup_size);
//             }
            
//             group_count[old_group_id] -= newgroup_size;
//             group_count[totgc] += newgroup_size;
            
//             totgc+=1;
            
//             ehh_before_norm += del_update;

//             //bool isDerivedGroup =  (!hm.hapData.hapEntries[locus].flipped && isDerived[ele.second[0]]) || (hm.hapData.hapEntries[locus].flipped && !isAncestral[ele.second[0]]); // just check first element to know if it is derived. 
//                 bool isDerivedGroup = isDerived[ele.second[0]];
//             if(isDerivedGroup) // if the core locus for this chr has 1 (derived), then update ehh1, otherwise ehh0
//             {
//                 ehh1_before_norm += del_update;
//             }else{
//                 ehh0_before_norm += del_update;
//             }
//         }
//         m.clear();
        

//         //printing 

//         double current_derived_ehh; 
//         double current_ancestral_ehh;
//         double current_ehh;
        
//         if(p.ALT){
//             current_ehh = ehh_before_norm*1.0/square_alt(n_c1+n_c0);
//             current_derived_ehh = ehh1_before_norm*1.0/square_alt(n_c1);
//             current_ancestral_ehh  = ehh0_before_norm*1.0/square_alt(n_c0);
            
//         }else{
//             current_ehh = ehh_before_norm*1.0/twice_num_pair(n_c1+n_c0);
//             current_derived_ehh = ehh1_before_norm*1.0/twice_num_pair(n_c1);
//             current_ancestral_ehh  = ehh0_before_norm*1.0/twice_num_pair(n_c0);
            
//         }

//         if(downstream){
//             (*fout) << std::fixed <<   -int(hm.mapData.mapEntries[locus].physicalPos -  hm.mapData.mapEntries[i].physicalPos)  << "\t"
//             <<  -(hm.mapData.mapEntries[locus].geneticPos -  hm.mapData.mapEntries[i].geneticPos)<< "\t"
//             << current_derived_ehh << "\t"
//             << current_ancestral_ehh << "\t";
//         }else{
//             (*fout) << std::fixed <<   hm.mapData.mapEntries[i].physicalPos -  hm.mapData.mapEntries[locus].physicalPos  << "\t"
//             <<  hm.mapData.mapEntries[i].geneticPos -  hm.mapData.mapEntries[locus].geneticPos<< "\t"
//             << current_derived_ehh << "\t"
//             << current_ancestral_ehh << "\t";
//         }

//         /*
//         cout<<"Iter "<<i-locus<<": EHH1["<<locus<<","<<i<<"]=";

//         for (int x = 0 ; x < totgc;  x++){
//             cout<<group_count[x]<<"("<<x<<")";
//         }

//         //print all elements of vector v
//         for (int x = 0 ; x < v.size();  x++){
//             cout<<v[x]<<" ";
//         }
        
//         cout<<endl;
//         */

            
//         if(unphased){
//             //(*fout) << current_notAncestral_ehh << "\t"
//             //        << current_notDerived_ehh << "\t";
//         }
//         (*fout) << current_ehh << endl;
//     }
// }


// /**
//  * @brief Calculate the EHH for a single locus
//  * @param query The query locus name
// */
// void EHH::calc_single_ehh(string query){
//     int numSnps = hm.mapData.nloci;
//     int numHaps = hm.hapData.nhaps;

//     //TODO
//     int locus = hm.queryFound(query);

//     unordered_map<int, vector<int> > m;
    
//     //calc_ehh_unidirection(locus, m, true); // downstream
//     calc_ehh_unidirection(locus, m, false); // upstream
// }

// // /**
// //  * @brief Calculate the EHH for a single locus
// //  * @param query The query locus phys pos
// // */
// // void EHH::calc_single_ehh(int query){
// //     int numSnps = hm.mapData.nloci;
// //     int numHaps = hm.hapData.nhaps;

// //     //TODO
// //     int locus = hm.queryFound(query);
    
// //     unordered_map<int, vector<int> > m;
// //     calc_ehh_unidirection(locus, m, false); // upstream
// //     calc_ehh_unidirection(locus, m, true); // downstream
    
// // }