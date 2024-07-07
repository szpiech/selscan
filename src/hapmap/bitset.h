#ifndef __BITSET_H__
#define __BITSET_H__

#include <cstdint>
#include <iostream>
#include <omp.h>

using namespace std;
class MyBitset{
    public:
    uint64_t* bits;
    int nbits;
    int nwords;
    int WORDSZ = 64;
    int num_1s = 0;

    MyBitset(int nbits){
        this->nbits = nbits;
        this->nwords = (nbits/WORDSZ) + 1; //idea: do ceil
        bits = new uint64_t[nwords];
        for (int i = 0; i < nwords; i++)
        {
            bits[i] = 0;
        }
    }

    // MyBitset(const MyBitset &source) {
    //     bits = new uint64_t;
    //     *bits = *source.bits;
    //     cout << "Copy constructor called" << endl;
    // }

    int count_1s(){
        int sum = 0;
        omp_set_num_threads(4);

        #pragma `omp parallel for reduction(+:sum)
        for (int k = 0; k < nwords; k++) {
            sum += __builtin_popcountll ((uint64_t)bits[k]);
        }
        return sum;
    }

    vector<int> get_bits(){
        uint64_t bitset;
        vector<int> bitvec;
        for (int k = 0; k < nwords; k++) {
            bitset = bits[k];
            //std::cout<<"B"<<bitset<<std::endl;
            
            while (bitset != 0) {
            uint64_t t = bitset & -bitset;
            int r = __builtin_ctzl(bitset);
            
            //std::cout<<(k * WORDSZ + r) << std::endl;
            
            
            bitvec.push_back(k * WORDSZ + r); //idea: reserve 1 counts
            //callback(k * WORDSZ + r);
            bitset ^= t;
            }
        }
        // for (int i = 0; i < nbits; i++)
        // {
        //     bitvec.push_back(get_bit(i));
        // }
        return bitvec;
    }

    void set_bit(int bit){
        bits[bit/WORDSZ] |= (uint64_t) 1 << (bit % WORDSZ);
    }

    void clear_bit(int bit){
        bits[bit/WORDSZ] &= ~((uint64_t) 1 << (bit % WORDSZ));
    }

    bool get_bit(int bit){
        return (bits[bit/WORDSZ] & ((uint64_t) 1 << (bit % WORDSZ))) != 0;
    }

    void print(){
        for (int i = 0; i < nbits; i++)
        {
            cout << get_bit(i);
        }
        cout << endl;
    }

    void print_pos(){
        for (int i = 0; i < nbits; i++)
        {
            if(get_bit(i)==1)
                cout << i << " ";
        }
        cout << endl;
    }

    MyBitset operator^(const MyBitset& b) {
        MyBitset xor_bitset(this->nbits);

        #pragma `omp parallel for
        for (int k = 0; k < this->nwords; k++) {
            xor_bitset.bits[k] = this->bits[k] ^ b.bits[k];
        }
        return xor_bitset;
    }

    ~MyBitset(){
        delete [] bits;
    }

};


#endif