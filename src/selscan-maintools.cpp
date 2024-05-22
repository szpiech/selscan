/* selscan -- a program to calculate EHH-based scans for positive selection in genomes
   Copyright (C) 2014-2024  Zachary A Szpiech
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
*/

#include "selscan-maintools.h"
pthread_mutex_t mutex_log = PTHREAD_MUTEX_INITIALIZER;

MainTools::MainTools(HapMap& hm, param_main& p,  ofstream* flog,  ofstream* fout){
    this->flog = flog;
    this->fout = fout; 
    this->hm = hm;
    this->p = p;
    //this->bar = new Bar();
    this->numThreads = p.numThreads;
}


/**
 * Calculate EHH in only one direction until cutoff is hit - upstream or downstream
*/
void EHH::calc_ehh_unidirection(int locus, unordered_map<unsigned int, vector<unsigned int> > & m, bool downstream){
    int numHaps = hm.hapData.nhaps;
    int numSnps = hm.mapData.nloci;
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
    vector<unsigned int> v = hm.hapData.hapEntries[locus].positions;
    
    //unordered_set<unsigned int> v = hm.all_positions[locus];
    // if(v.size()==0){
    //     n_c0 = numHaps;
    //     group_count[0] = numHaps;
    //     totgc+=1;
    //     ehh0_before_norm = twice_num_pair(n_c0);
    // }else if (v.size()==numHaps){ // all set
    //     group_count[0] = numHaps;
    //     totgc+=1;
    //     n_c1 = numHaps;
        
    //     for (int set_bit_pos : v){
    //         isDerived[set_bit_pos] = true;
    //     }
    //     ehh1_before_norm = twice_num_pair(n_c1);
    // }
    
    if(v.size() == 0 or  v.size() == numHaps){
        std::cerr<<"ERROR: Monomorphic site should not exist";
        throw 0;
    }

    hm.hapData.hapEntries[locus].flipped = false;
    if(hm.hapData.hapEntries[locus].flipped){
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
        if(ehh1_before_norm*1.0/twice_num_pair(n_c1) <= p.EHH_CUTOFF or ehh0_before_norm*1.0/twice_num_pair(n_c0)  <= p.EHH_CUTOFF){   // or cutoff, change for benchmarking against hapbin
            //std::cout<<"breaking"<<endl;
            break;
        }
        //*/
        double distance;
        
        if(downstream){
            distance = hm.mapData.mapEntries[i+1].physicalPos - hm.mapData.mapEntries[i].physicalPos;
        }else{
            distance = hm.mapData.mapEntries[i].physicalPos - hm.mapData.mapEntries[i-1].physicalPos;
        }

        // this should not happen as we already did integrity check previously
        if (distance < 0)
        {
            std::cerr << "ERROR: physical position not in ascending order.\n"; 
            throw 0;
        }
        
        
        // if(downstream){
        //     v = hm.hapData.hapEntries[i+1].xors;
        // }else{
        //     v = hm.hapData.hapEntries[i].xors;
        // }

        v = hm.hapData.hapEntries[i].positions;

        
        for (const unsigned int& set_bit_pos : v){
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

            //bool isDerivedGroup =  (!hm.hapData.hapEntries[locus].flipped && isDerived[ele.second[0]]) || (hm.hapData.hapEntries[locus].flipped && !isAncestral[ele.second[0]]); // just check first element to know if it is derived. 
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
            (*fout) << std::fixed <<   -int(hm.mapData.mapEntries[locus].physicalPos -  hm.mapData.mapEntries[i].physicalPos)  << "\t"
            <<  -(hm.mapData.mapEntries[locus].geneticPos -  hm.mapData.mapEntries[i].geneticPos)<< "\t"
            << current_derived_ehh << "\t"
            << current_ancestral_ehh << "\t";
        }else{
            (*fout) << std::fixed <<   hm.mapData.mapEntries[i].physicalPos -  hm.mapData.mapEntries[locus].physicalPos  << "\t"
            <<  hm.mapData.mapEntries[i].geneticPos -  hm.mapData.mapEntries[locus].geneticPos<< "\t"
            << current_derived_ehh << "\t"
            << current_ancestral_ehh << "\t";
        }

        cout<<"Iter "<<i-locus<<": EHH1["<<locus<<","<<i<<"]=";
            for (int x = 0 ; x < totgc;  x++){
                cout<<group_count[x]<<"("<<x<<")";
            }

            

            for (int x = 0 ; x < totgc;  x++){
                cout<<group_count[x]<<"("<<x<<")";
            }

            //print all elements of vector v
            for (int x = 0 ; x < v.size();  x++){
                cout<<v[x]<<" ";
            }
            
           cout<<endl;
        
        if(unphased){
            //(*fout) << current_notAncestral_ehh << "\t"
            //        << current_notDerived_ehh << "\t";
        }
        (*fout) << current_ehh << endl;
    }
}



string getOrder(uint64_t n_c2, uint64_t n_c1, uint64_t n_c0){
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


void IHS::updateEHH_from_split( map<int, vector<int> >& m, int* group_count, int* group_id, int& totgc, uint64_t ehh_before_norm[], uint64_t cehh_before_norm[], bool is1[], bool is2[]){
    for (const auto &ele : m) {
        int old_group_id = ele.first;
        int newgroup_size = ele.second.size() ;
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


void IHS::calc_ehh_unidirection_ihs_unphased(int locus, bool downstream){
    cout<<"exec unphased"<<endl;
    map<int, vector<int> > m1;
    map<int, vector<int> > m2;

    int numSnps = hm.hapData.nloci;
    int numHaps = hm.hapData.nhaps;

    uint64_t ehh_before_norm[3] = {0,0,0};
    uint64_t curr_ehh_before_norm[3] = {0,0,0};

    uint64_t cehh_before_norm[3] = {0,0,0};
    uint64_t curr_cehh_before_norm[3] = {0,0,0};


    uint64_t n_c[3] = {0,0,0};

    bool gap_skip = false;

    int group_count[numHaps];
    int group_id[numHaps];
    bool is0[numHaps];
    bool is1[numHaps];
    bool is2[numHaps];

    double* iHH[3] = {iHH0, iHH1, iHH2};
    double* ciHH[3] = {ciHH0, ciHH1, ciHH2};


    //will be vectorized with compile time flags
    for(int i = 0; i<numHaps; i++){
        group_count[i] = 0;
        group_id[i] = 0;
        is0[i] = false;
        is1[i] = false;
        is2[i] = false;
    }

    int totgc=0;
    vector<unsigned int> v1 = hm.hapData.hapEntries[locus].positions;
    vector<unsigned int> v2 = hm.hapData.hapEntries[locus].positions2;

    
    n_c[2] = v2.size();
    n_c[1] = v1.size();
    n_c[0] = numHaps - n_c[1] - n_c[2];
    string orderStr = getOrder(n_c[2], n_c[1], n_c[0]);
    
    

    for(int i = 0; i<3; i++){
        ehh_before_norm[i] = twice_num_pair(n_c[i]);
        group_count[i] = n_c[orderStr[i]-'0'];
    }
    
    //TODO
    if(group_count[0] == numHaps){ //monomorphic site
        totgc+=1;
    }else if(group_count[2] == 0 ){ //second clause is redundant
        if(group_count[1] + group_count[2] != numHaps){
            cerr<<"ERROR: gc2==0"<<endl;
            exit(1);
        }
        totgc+=2;
    }else{
        totgc+=3;
    }

    //group_count
    //[0] = most occurring
    //[1] = second most occurring
    //[2] = least occurring

    for (int set_bit_pos : v1){
        is1[set_bit_pos] = true;
        group_id[set_bit_pos] = 1;
    }

    for (int set_bit_pos : v2){
        is2[set_bit_pos] = true;
        group_id[set_bit_pos] = 2;
    }


    if(downstream){
        for(int i=0; i<3; i++){
            if(twice_num_pair(n_c[1])!=0){
                iHH[i][locus] += (curr_ehh_before_norm[i] + ehh_before_norm[i]) * 0.5 / twice_num_pair(n_c[i]);
            }
        }
    }
    
    for(int i=0; i<3; i++){
        curr_ehh_before_norm[i] = ehh_before_norm[i];
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
        
        
        //if(curr_ehh1_before_norm*1.0/n_c1_squared_minus < cutoff and curr_ehh0_before_norm*1.0/n_c0_squared_minus < cutoff){
        if(curr_ehh_before_norm[1]*1.0/twice_num_pair(n_c[1]) <= p.EHH_CUTOFF and curr_ehh_before_norm[0]*1.0/twice_num_pair(n_c[0])  <= p.EHH_CUTOFF){   // or cutoff, change for benchmarking against hapbin
            //std::cout<<"breaking"<<endl;
            break;
        }
        //TODO
        // if (current_derived_ehh <=  EHH_CUTOFF) dont add to ihh calculation, only do EHH ancestral 
        //     {
        //     }


        double distance;
        
        if(downstream){
            distance = hm.mapData.mapEntries[i+1].physicalPos - hm.mapData.mapEntries[i].physicalPos;
        }else{
            distance = hm.mapData.mapEntries[i].physicalPos - hm.mapData.mapEntries[i-1].physicalPos;
        }

        // this should not happen as we already did integrity check previously
        if (distance < 0)
        {
            std::cerr << "ERROR: physical position not in ascending order.\n"; 
            throw 0;
        }
        
        


        // if(distance> max_gap){
        //     gap_skip = true;
        //     break;
        // }
        if(distance > p.SCALE_PARAMETER){
            distance /= p.SCALE_PARAMETER;
        }
        //distance = 1; // for testing
        

        if(n_c[0] == numHaps or n_c[1] == numHaps or n_c[2] == numHaps){
            cerr<<"WARNING: Monomorphic site."<<endl;
        
            for(int i=0; i<3; i++){
                if(twice_num_pair(n_c[i])!=0){
                    iHH[i][locus] += (curr_ehh_before_norm[i] + ehh_before_norm[i]) * distance * 0.5 / twice_num_pair(n_c[i]);
                }
            }
            continue;
        }


        // ensure that in boundary we don't do any calculation
        // if(hm.hapData.hapEntries[i].positions.size() < v.size() && i!=numHaps-1 ){ 
        //     v = hm.hapData.hapEntries[i].positions;
        //     if(v.size()==0 or v.size()==numHaps){ // integrity check
        //         std::cerr<<"ERROR: Monomorphic site should not exist."<<endl;
        //         throw 0;
                
        //         if(twice_num_pair(n_c1)!=0){    // all 1s 
        //             iHH1[locus] += (curr_ehh1_before_norm + ehh1_before_norm) * distance * 0.5 / twice_num_pair(n_c1) ;
        //         }
        //         if(twice_num_pair(n_c0)!=0){   // all 0s 
        //             iHH0[locus] += (curr_ehh0_before_norm + ehh0_before_norm) * distance * 0.5 / twice_num_pair(n_c0)  ;
        //         }
        //         continue;
        //     }
        // }

        
        for (const unsigned int& set_bit_pos : v1){
            int old_group_id = group_id[set_bit_pos];
            m1[old_group_id].push_back(set_bit_pos);
        }

        for (const unsigned int& set_bit_pos : v2){
            int old_group_id = group_id[set_bit_pos];
            m2[old_group_id].push_back(set_bit_pos);
        }

        if(m1.size()==0 && m2.size()==0){
            //DO NOTHING
        }else if(m1.size()==0 || m2.size()==0){
            map<int, vector<int> >& m = (m1.size()==0) ? m2: m1;
            updateEHH_from_split(m, group_count, group_id, totgc, ehh_before_norm, cehh_before_norm, is1, is2);
        }else{ // both non-empty
            updateEHH_from_split(m1, group_count, group_id, totgc, ehh_before_norm,  cehh_before_norm,is1, is2);
            updateEHH_from_split(m2, group_count, group_id, totgc, ehh_before_norm,  cehh_before_norm,is1, is2);
        }

        for(int i=0; i<3; i++){
            if(twice_num_pair(n_c[i])!=0){
                if(ehh_before_norm[i]*1.0/twice_num_pair(n_c[i]) > p.EHH_CUTOFF){
                    iHH[i][locus] += (curr_ehh_before_norm[i] + ehh_before_norm[i]) * distance * 0.5 / twice_num_pair(n_c[i]);
                }
            }


            if(twice_num_pair(numHaps - n_c[i])!=0){
                if(cehh_before_norm[i]*1.0/twice_num_pair(numHaps - n_c[i]) > p.EHH_CUTOFF){
                    ciHH[i][locus] += (curr_cehh_before_norm[i] + cehh_before_norm[i]) * distance * 0.5 / twice_num_pair(numHaps - n_c[i]);
                }
            }

            curr_ehh_before_norm[i] = ehh_before_norm[i];
            curr_cehh_before_norm[i] = cehh_before_norm[i];


            //this shouldn't execute
            if(twice_num_pair(n_c[i])==0){
                curr_ehh_before_norm[i] = 0;
            }

            if(twice_num_pair(numHaps - n_c[i])==0){
                curr_cehh_before_norm[i] = 0;
            }
        }

        m1.clear();
        m2.clear();



        // if (nextLocus < 0 || nextLocus >= hm.mapData.nloci)
        // {
        //     pthread_mutex_lock(&mutex_log);
        //     (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << p.EHH_CUTOFF
        //             << ". ";
        //     if (!p.TRUNC){
        //         skipLocus = true;
        //         (*flog) << "Skipping calculation at position " << hm.mapData.mapEntries[locus].physicalPos << " id: " << hm.mapData.mapEntries[locus].locusName;
        //     }
        //     (*flog) << "\n";
        //     pthread_mutex_unlock(&mutex_log);
        //     break;
        // }
        // else if ( physicalDistance(nextLocus, downstream) > p.MAX_GAP )  
        // {
        //     pthread_mutex_lock(&mutex_log);
        //     (*flog) << "WARNING: Reached a gap of " << physicalDistance(nextLocus, downstream) << "bp > " << p.MAX_GAP 
        //     << "bp. Skipping calculation at position " <<  hm.mapData.mapEntries[locus].physicalPos << " id: " <<  hm.mapData.mapEntries[locus].locusName << "\n";
        //     pthread_mutex_unlock(&mutex_log);
        //     skipLocus = true;
        //     break;
        // }
    }

    
}

/**
 * Calculate EHH in only one direction until cutoff is hit - upstream or downstream
*/
void IHS::calc_ehh_unidirection_ihs(int locus, unordered_map<unsigned int, vector<unsigned int> > & m, bool downstream){
    int numSnps = hm.hapData.nloci;
    int numHaps = hm.hapData.nhaps;

    int total_iteration_of_m = 0;
    uint64_t ehh0_before_norm = 0;
    uint64_t ehh1_before_norm = 0;

    bool gap_skip = false;

    uint64_t curr_ehh0_before_norm = 0;
    uint64_t curr_ehh1_before_norm = 0;

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
    vector<unsigned int> v = hm.hapData.hapEntries[locus].positions;
    
    //unordered_set<unsigned int> v = hm.all_positions[locus];
    if(v.size()==0){
        n_c0 = numHaps;
        group_count[0] = numHaps;
        totgc+=1;
        ehh0_before_norm = twice_num_pair(n_c0);
    }else if (v.size()==numHaps){ // all set
        group_count[0] = numHaps;
        totgc+=1;
        n_c1 = numHaps;
        
        for (int set_bit_pos : v){
            isDerived[set_bit_pos] = true;
        }
        ehh1_before_norm = twice_num_pair(n_c1);
    }else{
        if(hm.hapData.hapEntries[locus].flipped){
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
    }

    if(downstream){
        // if(!calc_all)
        //     out_ehh<<"Iter "<<0<<": EHH1["<<locus<<","<<locus<<"]="<<1<<" "<<1<<endl;

        if(twice_num_pair(n_c1)!=0){
            iHH1[locus] += (curr_ehh1_before_norm + ehh1_before_norm) * 0.5 / twice_num_pair(n_c1);
            //cout<<"Summing "<<1.0*curr_ehh1_before_norm/n_c1_squared_minus<<" "<<1.0*ehh1_before_norm/n_c1_squared_minus<<endl;
        }
        if(twice_num_pair(n_c0)!=0){
            iHH0[locus] += (curr_ehh0_before_norm + ehh0_before_norm) * 0.5 / twice_num_pair(n_c0);
        }
    }
    
    curr_ehh1_before_norm = ehh1_before_norm;
    curr_ehh0_before_norm = ehh0_before_norm;

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
        if(curr_ehh1_before_norm*1.0/twice_num_pair(n_c1) <= p.EHH_CUTOFF and curr_ehh0_before_norm*1.0/twice_num_pair(n_c0)  <= p.EHH_CUTOFF){   // or cutoff, change for benchmarking against hapbin
            //std::cout<<"breaking"<<endl;
            break;
        }
        //TODO
        // if (current_derived_ehh <=  EHH_CUTOFF) dont add to ihh calculation, only do EHH ancestral 
        //     {
        //     }


        double distance;
        
        if(downstream){
            distance = hm.mapData.mapEntries[i+1].physicalPos - hm.mapData.mapEntries[i].physicalPos;
        }else{
            distance = hm.mapData.mapEntries[i].physicalPos - hm.mapData.mapEntries[i-1].physicalPos;
        }

        // this should not happen as we already did integrity check previously
        if (distance < 0)
        {
            std::cerr << "ERROR: physical position not in ascending order.\n"; 
            throw 0;
        }
        
        


        // if(distance> max_gap){
        //     gap_skip = true;
        //     break;
        // }
        if(distance > p.SCALE_PARAMETER){
            distance /= p.SCALE_PARAMETER;
        }
        //distance = 1; // for testing
        
        if(downstream){
            v = hm.hapData.hapEntries[i+1].xors;
        }else{
            v = hm.hapData.hapEntries[i].xors;
        }

        // ensure that in boundary we don't do any calculation
        if(hm.hapData.hapEntries[i].positions.size() < v.size() && i!=numHaps-1 ){ 
            v = hm.hapData.hapEntries[i].positions;
            if(v.size()==0 or v.size()==numHaps){ // integrity check
                std::cerr<<"ERROR: Monomorphic site should not exist."<<endl;
                throw 0;
                
                if(twice_num_pair(n_c1)!=0){    // all 1s 
                    iHH1[locus] += (curr_ehh1_before_norm + ehh1_before_norm) * distance * 0.5 / twice_num_pair(n_c1) ;
                }
                if(twice_num_pair(n_c0)!=0){   // all 0s 
                    iHH0[locus] += (curr_ehh0_before_norm + ehh0_before_norm) * distance * 0.5 / twice_num_pair(n_c0)  ;
                }
                continue;
            }
        }

        v = hm.hapData.hapEntries[i].positions;
        
        for (const unsigned int& set_bit_pos : v){
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
            
            bool isDerivedGroup =  (!hm.hapData.hapEntries[locus].flipped && isDerived[ele.second[0]]) || (hm.hapData.hapEntries[locus].flipped && !isAncestral[ele.second[0]]); // just check first element to know if it is derived. 
                
            if(isDerivedGroup) // if the core locus for this chr has 1 (derived), then update ehh1, otherwise ehh0
            {
                ehh1_before_norm += del_update;
            }else{
                ehh0_before_norm += del_update;
            }
        }

        if(twice_num_pair(n_c1)!=0){

            if(ehh1_before_norm*1.0/twice_num_pair(n_c1) > p.EHH_CUTOFF){
               
                iHH1[locus] += (curr_ehh1_before_norm + ehh1_before_norm) * distance * 0.5 / twice_num_pair(n_c1);
            }
            
            //cout<<"Summing "<<1.0*curr_ehh1_before_norm/n_c1_squared_minus<<" "<<1.0*ehh1_before_norm/n_c1_squared_minus<<endl;
        }

        if(twice_num_pair(n_c0)!=0){
            if(ehh0_before_norm*1.0/twice_num_pair(n_c0) > p.EHH_CUTOFF){
                iHH0[locus] += (curr_ehh0_before_norm + ehh0_before_norm) * distance * 0.5 / twice_num_pair(n_c0);

            }
        }

        curr_ehh1_before_norm = ehh1_before_norm;
        curr_ehh0_before_norm = ehh0_before_norm;

        // ehh1[locus] = 1.0*ehh1_before_norm/n_c1_squared_minus;
        // ehh0[locus] = 1.0*ehh0_before_norm/n_c0_squared_minus;

        //this shouldn't execute
        if(twice_num_pair(n_c1)==0){
            curr_ehh1_before_norm = 0;
        }
        if(twice_num_pair(n_c0)==0){
            curr_ehh0_before_norm = 0;
        }

        m.clear();



        // if (nextLocus < 0 || nextLocus >= hm.mapData.nloci)
        // {
        //     pthread_mutex_lock(&mutex_log);
        //     (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << p.EHH_CUTOFF
        //             << ". ";
        //     if (!p.TRUNC){
        //         skipLocus = true;
        //         (*flog) << "Skipping calculation at position " << hm.mapData.mapEntries[locus].physicalPos << " id: " << hm.mapData.mapEntries[locus].locusName;
        //     }
        //     (*flog) << "\n";
        //     pthread_mutex_unlock(&mutex_log);
        //     break;
        // }
        // else if ( physicalDistance(nextLocus, downstream) > p.MAX_GAP )  
        // {
        //     pthread_mutex_lock(&mutex_log);
        //     (*flog) << "WARNING: Reached a gap of " << physicalDistance(nextLocus, downstream) << "bp > " << p.MAX_GAP 
        //     << "bp. Skipping calculation at position " <<  hm.mapData.mapEntries[locus].physicalPos << " id: " <<  hm.mapData.mapEntries[locus].locusName << "\n";
        //     pthread_mutex_unlock(&mutex_log);
        //     skipLocus = true;
        //     break;
        // }
    }

    
}



void IHS::thread_ihs(int tid, unordered_map<unsigned int, vector<unsigned int> >& m, unordered_map<unsigned int, vector<unsigned int> >& md, IHS* ehh_obj){
    int elem_per_block = floor(ehh_obj->hm.hapData.nloci/ehh_obj->numThreads);
    int start = tid*elem_per_block ;
    int end = start + elem_per_block  ;
    if(tid == ehh_obj->numThreads-1 ){
        end = ehh_obj->hm.hapData.nloci;
    }

    //#pragma omp parallel 
    for(int locus = start; locus< end; locus++){
        ehh_obj->calc_ihh(locus);
    }
    
    //ehh_obj->log_string_per_thread[tid]+="finishing thread #"+to_string(tid)+"\n"; 
    
    //Logger::write("finishing thread # "+to_string(tid)+" at "+to_string(readTimer())+"\n");
    pthread_mutex_lock(&mutex_log);
    (*(ehh_obj->flog))<<("finishing thread # "+to_string(tid)+" at "+to_string(MainTools::readTimer())+"\n");
    pthread_mutex_unlock(&mutex_log);

}

void IHS::ihs_main(){
    iHH0 = new double[hm.hapData.nloci];
    iHH1 = new double[hm.hapData.nloci];

    if(hm.hapData.unphased){
        iHH2 = new double[hm.hapData.nloci];
        ciHH0 = new double[hm.hapData.nloci];
        ciHH1 = new double[hm.hapData.nloci];
        ciHH2 = new double[hm.hapData.nloci];
    }
    
    std::unordered_map<unsigned int, std::vector<unsigned int> > map_per_thread[numThreads];
    std::unordered_map<unsigned int, std::vector<unsigned int> > mapd_per_thread[numThreads];


    bool openmp_enabled = false;
    // two different ways to parallelize: first block does pthread, second block does openmp
    if (!openmp_enabled)
    {
        int total_calc_to_be_done = hm.hapData.nloci;
        cout<<"total calc to be done: "<<total_calc_to_be_done<<endl;
        std::thread *myThreads = new std::thread[numThreads];
        for (int i = 0; i < numThreads; i++)
        {
            myThreads[i] = std::thread(thread_ihs, i, std::ref(map_per_thread[i]),  std::ref(mapd_per_thread[i]), this);
        }
        for (int i = 0; i < numThreads; i++)
        {
            myThreads[i].join(); // Join will block our main thread, and so the program won't exit until all finish
        }
        delete[] myThreads;


        //Logger::write("all threads finished. now calculating ihh...\n");
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
    delete[] iHH0;
    delete[] iHH1; 
    if(hm.hapData.unphased){
        delete[] iHH2;
        delete[] ciHH0;
        delete[] ciHH1;
        delete[] ciHH2;
    }

    return;
}


/**
 * @brief Calculate unphased IHS
 * @param locus The locus 
*/
double IHS::calc_ihs_unphased(int locus){
     double iHS2 = log10(iHH2[locus]/ciHH2[locus]);
     double iHS0 = log10(iHH0[locus]/ciHH0[locus]);
    if(iHS2 > iHS0){
        return iHS2;
    }else{
        return -iHS0;
    }
}



/**
 * @brief Calculate the EHH for a single locus as part of IHH routine
 * @param locus The locus 
*/
void IHS::calc_ihh(int locus){
    int numSnps = hm.mapData.nloci;
    int numHaps = hm.hapData.nhaps;
    
    unordered_map<unsigned int, vector<unsigned int> > m;
    iHH0[locus] = 0;
    iHH1[locus] = 0;

    if(hm.hapData.unphased){
        calc_ehh_unidirection_ihs_unphased(locus, false); // upstream
        calc_ehh_unidirection_ihs_unphased(locus, true); // downstream
    }else{
        calc_ehh_unidirection_ihs(locus, m, false); // upstream
        calc_ehh_unidirection_ihs(locus, m, true); // downstream
    }
    
    if(!hm.hapData.unphased){
        //handle all corner cases
        if(hm.hapData.hapEntries[locus].positions.size()==0){
            iHH1[locus] = 1;
        }
        if(hm.hapData.hapEntries[locus].positions.size()==numHaps){ //iHH0[locus]==0
            iHH0[locus] = 1;
        }
        if(hm.hapData.hapEntries[locus].positions.size()==1){
            if(locus == 0 or locus == numSnps-1){
                iHH1[locus] = 0.5;
            }else{
                iHH1[locus] = 1;
            }
        }
        if(hm.hapData.hapEntries[locus].positions.size()==numHaps-1){
            if(locus == 0 or locus == numSnps-1){
                iHH0[locus] = 0.5;
            }else{
                iHH0[locus] = 1;
            }
        }
        //PRINTHERE
        //  if(hm.getMAF(i) >= min_maf && 1-hm.getMAF(i) >= min_maf){
        //         sprintf(str, "%d %d %f %f %f %f\n", hm.mentries[i].phyPos, hm.mentries[i].locId, hm.all_positions[i].size()*1.0/numHaps, iHH1[i], iHH0[i], log10(iHH1[i]/ iHH0[i]));
        //         out_ihs<<str;
        //      }

        (*fout) << std::fixed <<   hm.mapData.mapEntries[locus].locusName << " " <<   hm.mapData.mapEntries[locus].physicalPos << " "
                <<  hm.mapData.mapEntries[locus].locId << " " << hm.hapData.calcFreq(locus) << " "
                << iHH1[locus] << " " << iHH0[locus] <<" "<< iHH1[locus]*1.0/iHH0[locus] <<endl;
        
            
    }

    if(hm.hapData.unphased){
        (*fout) << std::fixed <<   hm.mapData.mapEntries[locus].locusName << " " <<   hm.mapData.mapEntries[locus].physicalPos << " "
                <<  hm.mapData.mapEntries[locus].locId << " " << hm.hapData.calcFreq(locus) << " "
                << iHH1[locus] << " " << iHH0[locus] <<" "<< calc_ihs_unphased(locus) <<endl;
    }
    
}


/**
 * @brief Calculate the EHH for a single locus
 * @param query The query locus name
*/
void EHH::calc_single_ehh(string query){
    int numSnps = hm.mapData.nloci;
    int numHaps = hm.hapData.nhaps;

    //TODO
    int locus = hm.queryFound(query);
    cout<<"LOP "<<locus<<endl;
    unordered_map<unsigned int, vector<unsigned int> > m;
    
    //calc_ehh_unidirection(locus, m, true); // downstream
    calc_ehh_unidirection(locus, m, false); // upstream
    
    
}

// /**
//  * @brief Calculate the EHH for a single locus
//  * @param query The query locus phys pos
// */
// void EHH::calc_single_ehh(unsigned int query){
//     int numSnps = hm.mapData.nloci;
//     int numHaps = hm.hapData.nhaps;

//     //TODO
//     int locus = hm.queryFound(query);
    
//     unordered_map<unsigned int, vector<unsigned int> > m;
//     calc_ehh_unidirection(locus, m, false); // upstream
//     calc_ehh_unidirection(locus, m, true); // downstream
    
// }

void XPIHH::xpihh_main()
{
    int nloci = hm.mapData.nloci;
    ihh_p1 = new double[nloci];
    ihh_p2 = new double[nloci];

    barInit(*bar, hm.mapData.nloci, 78); //REWRITE

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
    
    hm.hapData.releaseHapData();
    hm.hapData2.releaseHapData();
    
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
    int numSnps = obj->hm.hapData.nhaps;
    int elem_per_block = floor(numSnps/obj->numThreads);
    int start = tid*elem_per_block ;
    int end = start + elem_per_block  ;
    if(tid == obj->numThreads-1 ){
        end = numSnps;
    }

    int step = (numSnps  / obj->numThreads) / (obj->bar->totalTicks);
    if (step == 0) step = 1;
    //if total 20 tasks: and 4 threads: t0: 0, 4, 8, 12, 16: t1: 1, 5, 9, 13, 17  
    //for (int locus = tid; locus < hm.mapData.nloci; locus += numThreads)
    //#pragma omp parallel 
    for(int locus = start; locus< end; locus++){
        if (locus % step == 0) advanceBar(*(obj->bar), double(step));
        obj->calc_xpihh(locus);
    }

    pthread_mutex_lock(&mutex_log);
    (*(obj->flog))<<("finishing thread # "+to_string(tid)+" at "+to_string(MainTools::readTimer())+"\n");
    pthread_mutex_unlock(&mutex_log);
}


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

    vector<unsigned int> twos_p1 = hm.hapData.hapEntries[locus].positions2;
    vector<unsigned int> twos_p2 = hm.hapData2.hapEntries[locus].positions2;

    
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

