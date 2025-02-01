#ifndef __XPIHH_H__
#define __XPIHH_H__

#include "selscan-stats.h"
#include <unordered_map>

using namespace std;
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

class XPIHH_ehh_data{
    public:
    double prev_ehh_before_norm = 0;
    double curr_ehh_before_norm = 0;
    int n_c[2] = {0,0};

    int nhaps;
    int totgc = 0;
    int* group_count;
    // bool* is1;
    int* group_id;
    double normalizer;

    ~XPIHH_ehh_data();
    // inline int square_alt(int n);
    // inline double twice_num_pair(int n);
    // void init(int nhaps, const vector <int>* positions, bool ALT, bool WAGH );
    // void init_for_pooled(const vector <int>* positions,const vector <int>* positions2, int nhaps_p1, bool ALT, bool WAGH );
    // void initialize_core(const vector <int>* v);
    // void initialize_core_pooled(const vector <int>* v, const vector <int>* v2,  int nhaps_p1); 
}; 

class XPIHH_ehh_data_unphased: public XPIHH_ehh_data{
    public:
    int n_c[3] = {0,0,0}; //  overriden

    string getOrder(int n_c2, int n_c1, int n_c0){
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
    void assign_groups(MyBitset* all1s, MyBitset* all2s, MyBitset* second_all1s = NULL, MyBitset* second_all2s = NULL, int nhaps_first = -1){

        string orderStr = getOrder(n_c[2], n_c[1], n_c[0]);
        for(int i = 0; i<3; i++){
            // curr_ehh_before_norm[i] = normalizer[i];
            // curr_cehh_before_norm[i] = normalizer_not[i];
            group_count[i] = n_c[orderStr[i]-'0']; // most occurring gets 0 id, 
            //group_core[i] = orderStr[i]-'0';
        }

        //group_count
        //[0] = most occurring
        //[1] = second most occurring
        //[2] = least occurring

        int pos_of_012[3] = {0,0,0};
        for(int i = 0; i<3; i++){
            pos_of_012[orderStr[i]-'0'] = i ;
        }

        for (int i = 0; i<nhaps; i++){
            group_id[i] = pos_of_012[0];
        }

        ACTION_ON_ALL_SET_BITS(all1s, {
            //is1[set_bit_pos] = true;
            group_id[set_bit_pos] = pos_of_012[1];   
        });

        ACTION_ON_ALL_SET_BITS(all2s, {
            //is2[set_bit_pos] = true;
            group_id[set_bit_pos] = pos_of_012[2]; 
        });


        if(nhaps_first != -1){ // pooled
            ACTION_ON_ALL_SET_BITS(second_all1s, {
                group_id[set_bit_pos  + nhaps_first] = pos_of_012[1];   
            });

            ACTION_ON_ALL_SET_BITS(second_all2s, {
                group_id[set_bit_pos + nhaps_first] = pos_of_012[2]; 
            });
        }

        if(group_count[0] == nhaps){ //monomorphic site
            totgc+=1;
        }else if(group_count[2] == 0 ){ //second clause is redundant
            if(group_count[0] + group_count[1] != nhaps){
                cerr<<"ERROR: gc2==0"<<endl;
                exit(2);
            }
            totgc+=2;
        }else{
            totgc+=3;
        }

        


    //     double normalizer[3];
    // double normalizer_not[3];
    //  for(int i = 0; i<3; i++){
    //     normalizer[i] = twice_num_pair_or_square(n_c[i], p  .ALT);
    //     normalizer_not[i] = twice_num_pair_or_square(numHaps - n_c[i], p.ALT); 
    //  }
    }



    // ~XPIHH_ehh_data_unphased();
    // void init(int nhaps, const vector <int>* positions, bool ALT, bool WAGH );
    // void init_for_pooled(const vector <int>* positions,const vector <int>* positions2, int nhaps_p1, bool ALT, bool WAGH );
    // void initialize_core(const vector <int>* v);
    // void initialize_core_pooled(const vector <int>* v, const vector <int>* v2,  int nhaps_p1); 
};

class XPIHH : public SelscanStats{
    public:
        XPIHH(const std::unique_ptr<HapMap>&  hm, param_main& params) : SelscanStats(hm, params){}
        void main();
        ~XPIHH(){}

    private:
        static pthread_mutex_t mutex_log;
        int max_extend;
        pair<double, double> calc_xpihh(int locus);
        pair<double, double> calc_ehh_unidirection(int locus, bool downstream);
        pair<double, double> calc_ehh_unidirection_unphased(int locus, bool downstream);

        void updateEHH_from_split(const unordered_map<int, vector<int> > & m, XPIHH_ehh_data* p);
        void updateGroup_from_split_unphased( unordered_map<int, vector<int> >& m, int* group_count, int* group_id, int& totgc);
};


#endif