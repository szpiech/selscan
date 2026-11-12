#include <memory>
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
    //int prev_index = locus;
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
            output[i].pdist = downstream? -physicalDistance(i,locus, downstream): physicalDistance(i,locus, downstream);
            output[i].gdist = downstream? -geneticDistance(i,locus, downstream): geneticDistance(i,locus, downstream);
            output[i].print = true;
            //cout<<i<< " "<<output[i].pdist<<" "<<output[i].ehh0<<" "<<output[i].ehh1<<" "<<output[i].ehh<<endl;
        }else{
            // PHASE 2: GET MAP
            GroupSplitMap& m = threadSplitMap(0);   // P5
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
            output[i].pdist = downstream? -physicalDistance(i,locus, downstream): physicalDistance(i,locus, downstream);
            output[i].gdist = downstream? -geneticDistance(i,locus, downstream): geneticDistance(i,locus, downstream);
            output[i].print = true;
            //cout<<i<< " "<<output[i].pdist<<" "<<output[i].ehh0<<" "<<output[i].ehh1<<" "<<output[i].ehh<<endl;

        }

        // if(totgc == numHaps) {
        //     //std::cerr<<"Break reason for locus "<<locus<<":: ALL_UNIQUE."<<endl;
        //     break;
        // }

        if(physicalDistance(i,locus, downstream) > p.QWIN) {
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
void EHH::main(string query_list){
    
    init_output_and_querymap();
    //init_global_fout("ehh",query);

    int numSnps = hm->mapData->nloci;

    std::vector<std::string> tokens;
    std::stringstream ss(query_list);
    std::string token;

    while (std::getline(ss, token, ',')) {
        if (token.empty()) {
            throw std::invalid_argument("Empty token found in query string.");
        }

        tokens.push_back(token);
    }

    // Extra check: no trailing comma
    if (!query_list.empty() && query_list.back() == ',') {
        throw std::invalid_argument("Query string ends with a comma.");
    }

    
    for(const string& query: tokens){

        int locus = this->queryFound(query);
        if(locus == -1){
            continue; // skip if query locus not found
        }

        string outFilename = get_filename_base("ehh", query);
        ofstream ehhOUT(outFilename+".out");

        if(p.UNPHASED){
            for (int i = 0; i < numSnps; i++){
                output[i].print = false;   // rows of a previous locus (B30)
            }
            calc_ehh_unidirection_unphased(locus, true); // downstream
            calc_ehh_unidirection_unphased(locus, false); // upstream

            ehhOUT << std::fixed <<   "chr\tpdist\tgdist\tehh2\tehh0\tcehh2\tcehh0\tehh" << "\n";
            for (int i = 0; i < numSnps; i++){
                if(!output[i].print){
                    continue;
                }
                ehhOUT << std::fixed <<   hm->mapData->mapEntries[i].chr << "\t"
                <<   output[i].pdist  << "\t"
                <<  output[i].gdist << "\t"
                << output[i].ehh2 << "\t"
                << output[i].ehh0  << "\t"
                << output[i].cehh2  << "\t"
                << output[i].cehh0  << "\t"
                << output[i].ehh  << "";
                ehhOUT << endl;
            }
        
        }else{
            for (int i = 0; i < numSnps; i++){
                (output[i].print = false);
            }

            calc_ehh_unidirection(locus, true); // downstream
            calc_ehh_unidirection(locus, false); // upstream

            ehhOUT << std::fixed <<   "chr\tpdist\tgdist\tderEHH\tancEHH\tEHH" << "\n";
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
                ehhOUT << std::fixed <<   hm->mapData->mapEntries[i].chr << "\t"
                <<   output[i].pdist  << "\t"
                <<  output[i].gdist << "\t"
                << output[i].ehh1 << "\t"
                << output[i].ehh0  << "\t"
                << output[i].ehh  << "";
                ehhOUT << endl;
            }
        }
        ehhOUT.close();
    }

    




    
}


void EHH::calc_ehh_unidirection_unphased(int locus, bool downstream){
    // B30: the table written by --ehh --unphased. At every site it shows the
    // current EHH of the two homozygous classes (ehh2, ehh0), their
    // complements (cehh2: individuals not homozygous 2; cehh0: not
    // homozygous 0), and whole-sample EHH. 3.0.1 never set the complement,
    // printed 0 for whole-sample EHH at the core, froze ehh2/ehh0 once they
    // fell below --cutoff (that gate belongs to the iHS integration, not to
    // the decay table), printed nan for a class of fewer than 2, and leaked
    // its work arrays on the early returns.

    GroupSplitMap& m = threadSplitMap(0);   // P5

    const int numHaps = hm->hapData->nhaps;

    std::unique_ptr<int[]> group_count_own(new int[numHaps]()), group_id_own(new int[numHaps]()), group_core_own(new int[numHaps]());
    std::unique_ptr<bool[]> is1_own(new bool[numHaps]()), is2_own(new bool[numHaps]());
    int* group_count = group_count_own.get();
    int* group_id = group_id_own.get();
    int* group_core = group_core_own.get();
    bool* is1 = is1_own.get();
    bool* is2 = is2_own.get();

    int n_c[3] = {0,0,0};
    n_c[1] = hm->hapData->get_n_c1(locus);
    n_c[2] = hm->hapData->get_n_c2(locus);
    n_c[0] = numHaps - n_c[1] - n_c[2]; // assume no missing
    if(n_c[1] + n_c[2] + n_c[0] != numHaps){
        HANDLE_ERROR("n_c1 + n_c2 + n_c0 != numHaps");
    }

    string orderStr = getOrder(n_c[2], n_c[1], n_c[0]);

    double normalizer[3];
    double normalizer_not[3];
    const double normalizer_just = twice_num_pair_or_square(numHaps, p.ALT);
    for(int i = 0; i<3; i++){
        normalizer[i] = twice_num_pair_or_square(n_c[i], p.ALT);
        normalizer_not[i] = twice_num_pair_or_square(numHaps - n_c[i], p.ALT);
    }
    auto ratio = [](double num, double den){ return den == 0 ? double(UNDEFINED_EHH) : num / den; };

    // group_count[0..2]: most, second most and least frequent genotype class
    for(int i = 0; i<3; i++){
        group_count[i] = n_c[orderStr[i]-'0'];
        group_core[i] = orderStr[i]-'0';
    }
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

    int totgc = 0;
    if(group_count[0] == numHaps){ //monomorphic site
        totgc+=1;
    }else if(group_count[2] == 0 ){
        if(group_count[0] + group_count[1] != numHaps){
            HANDLE_ERROR("gc2==0 locus: "+to_string(locus));
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
    if (freqHetGT > 1-p.MAF)
    {
        (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
                << " (number " << locus << ") has too many hets. Skipping calculation at this locus. "
                << "het: " <<  n_c[1] << " hom0: " << n_c[0] << " hom1: " << n_c[2] << ".\n";
        return;
    }
    if (normalizer[2] == 0 || normalizer[0] == 0) {
        (*flog) << "WARNING: locus " << hm->mapData->mapEntries[locus].locusName
                << " has a homozygous class of fewer than 2 individuals. Undefined EHH is reported as " << UNDEFINED_EHH << ".\n";
    }

    // core: each class and its complement start at 1, as in the unphased
    // iHS kernel; whole-sample EHH is the homozygosity of the 3 classes
    output[locus].ehh2 = ratio(normalizer[2], normalizer[2]);
    output[locus].ehh0 = ratio(normalizer[0], normalizer[0]);
    output[locus].cehh2 = ratio(normalizer_not[2], normalizer_not[2]);
    output[locus].cehh0 = ratio(normalizer_not[0], normalizer_not[0]);
    output[locus].ehh = ratio(normalizer[0] + normalizer[1] + normalizer[2], normalizer_just);
    output[locus].gdist = 0;
    output[locus].pdist = 0;
    output[locus].print = true;

    int i = locus;  // locus == core_locus
    while(true){
        if(nextLocOutOfBounds(i, downstream)) break;
        i = (downstream) ? i-1 : i+1;

        ACTION_ON_ALL_SET_BITS(get_all_2s(i), {
            int old_group_id = group_id[set_bit_pos];
            m[old_group_id].push_back(set_bit_pos);
        });
        updateEHH_from_split_unphased(m, group_count, group_id, totgc, is1, is2, group_core);
        m.clear();

        ACTION_ON_ALL_SET_BITS(get_all_1s(i), {
            int old_group_id = group_id[set_bit_pos];
            m[old_group_id].push_back(set_bit_pos);
        });
        updateEHH_from_split_unphased(m, group_count, group_id, totgc, is1, is2, group_core);
        m.clear();

        // homozygosity of each class, of each complement and of the sample
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

        output[i].ehh2 = ratio(ehh2_before_norm, normalizer[2]);
        output[i].ehh0 = ratio(ehh0_before_norm, normalizer[0]);
        output[i].cehh2 = ratio(cehh2_before_norm, normalizer_not[2]);
        output[i].cehh0 = ratio(cehh0_before_norm, normalizer_not[0]);
        output[i].ehh = ratio(ehh_before_norm, normalizer_just);
        output[i].pdist = downstream? -physicalDistance(i,locus, downstream): physicalDistance(i,locus, downstream);
        output[i].gdist = downstream? -geneticDistance(i,locus, downstream): geneticDistance(i,locus, downstream);
        output[i].print = true;

        // same window rule as the phased table
        if(physicalDistance(i,locus, downstream) > p.QWIN) {
            output[i].print = false;
            break;
        }
    }
}
