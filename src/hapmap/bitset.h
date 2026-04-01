#ifndef __BITSET_H__
#define __BITSET_H__

#include <cstdint>
#include <iostream>
// #include <omp.h>
#include <vector>
#include <cstdlib>   // For std::aligned_alloc
#include <cstring>   // For std::memset
#include <fstream>

#ifdef _WIN32
    #include <intrin.h> // For __popcnt64 and _tzcnt_u64
// #else
//     #include <x86intrin.h> // For __builtin_popcountll and __builtin_ctzl
#endif


using namespace std;


class MyBitset{
    public:
    uint64_t* bits;
    int nbits;  // number of bits
    int nwords;
    static constexpr int WORDSZ = 64;
    int num_1s = 0;

    MyBitset(int nbits){
        this->nbits = nbits;
        this->nwords = (nbits + WORDSZ - 1) / WORDSZ; // true ceil
        bits = (uint64_t*)calloc(nwords, sizeof(uint64_t));
        if (!bits) {
            cerr << "ERROR: Memory allocation failed for MyBitset with nbits = " << nbits << endl;
            exit(EXIT_FAILURE);
        }
    }

    //untested
    /*
    MyBitset(const MyBitset& other) {
        //     cout << "Copy constructor called" << endl;
        nbits = other.nbits;
        nwords = other.nwords;
        bits = new uint64_t[nwords];
        for (int i = 0; i < nwords; ++i) {
            bits[i] = other.bits[i];
        }
    }

    // Assignment operator
    MyBitset& operator=(const MyBitset& other) {
        if (this == &other) {
            return *this;
        }
        delete[] bits;

        nbits = other.nbits;
        nwords = other.nwords;

        bits = new uint64_t[nwords];
        for (int i = 0; i < nwords; ++i) {
            bits[i] = other.bits[i];
        }
        return *this;
    }
    // XOR operator
    MyBitset operator^(const MyBitset& other) const {
        if (nbits != other.nbits) {
            throw std::invalid_argument("Bitsets must be of the same size for XOR operation");
        }
        MyBitset result(nbits);
        for (int i = 0; i < nbits; ++i) {
            result.bits[i] = bits[i] ^ other.bits[i];
        }
        return result;
    }
    // MyBitset operator^(const MyBitset& b) { //wrong?
    //     MyBitset xor_bitset(this->nbits);

    //     //#pragma `omp parallel for
    //     for (int k = 0; k < this->nwords; k++) {
    //         xor_bitset.bits[k] = this->bits[k] ^ b.bits[k];
    //     }
    //     return xor_bitset;
    // }
    */

    int count_1s(){
        int sum = 0;
        //omp_set_num_threads(4);

        //#pragma `omp parallel for reduction(+:sum)
        for (int k = 0; k < nwords; k++) {
            sum += __builtin_popcountll ((uint64_t)bits[k]);
        }

        
        return sum;
    }

    int count_1s_fast(){
        int sum = 0;
        uint64_t* p = bits;
        uint64_t* end = bits + nwords;

        while (p != end)
            sum += __builtin_popcountll(*p++);

        return sum;
    }


    // // Declare a function pointer type
    // typedef int (*ActionSetBitPtr)(int);
    // void static performOperation(ActionSetBitPtr opFunc, int x) {
    //     opFunc(x);
    // }

    // Perform an action on all set bits
    /*
    template <typename Action>
    void performActionOnSetBits(Action action) {
        for (int k = 0; k < nwords; k++) {
            uint64_t bitset = bits[k];
            while (bitset != 0) {
                uint64_t t = bitset & -bitset; // Get the lowest set bit
                int r = __builtin_ctzl(bitset); // Get the index of the lowest set      
                int set_bit_pos = (k * WORDSZ + r); // Calculate the position of the set bit
                bitset ^= t; // Clear the lowest set bit
                action(set_bit_pos); // Perform the action on the set bit position
            }
        }
    }
    void action_on_all_set_bits(){
        for (int k = 0; k < nwords; k++) {
            uint64_t bitset = bits[k];
            while (bitset != 0) {
                uint64_t t = bitset & -bitset;
                int r = __builtin_ctzl(bitset);
                // DO ACTION HERE
                bitset ^= t;
            }
        }
    }
    */

    

    void set_bits(vector <unsigned int>& v){
        uint64_t bitset;
        v.resize(num_1s);
        int i = 0;
        for (int k = 0; k < nwords; k++) {
            bitset = bits[k];
            while (bitset != 0) {
                uint64_t t = bitset & -bitset;
                int r = __builtin_ctzl(bitset);
                v[i++] = k * WORDSZ + r; //idea: reserve 1 counts
                bitset ^= t;
            }
        }
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
        if(bit > nbits-1){
            cerr << "Error: bit out of range: " << bit<< " "<<nbits-1 << endl;
            throw (EXIT_FAILURE);
        }
        //bits[bit/WORDSZ] |= (uint64_t) 1 << (bit % WORDSZ);
        bits[bit >> 6] |= (uint64_t)1 << (bit & 63);
    }

    void clear_bit(int bit){
        bits[bit/WORDSZ] &= ~((uint64_t) 1 << (bit % WORDSZ));
    }

    bool get_bit(int bit){
        if(bit > nbits-1){
            cerr << "Error: bit out of range: " << bit<< " "<<nbits-1 << endl;
            throw (EXIT_FAILURE);
        }
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

    void print_pos_to_file(ofstream* fout){
        for (int i = 0; i < nbits; i++)
        {
            if(get_bit(i)==1)
                (*fout) << i << " ";
        }
        (*fout) << endl;
    }

    // Copy only the raw bit data into an already-initialized bitset.
    // Assumes this bitset already has the correct size.
    void copy_bits_from(const MyBitset& other){
        if (nbits != other.nbits || nwords != other.nwords) {
            cerr << "ERROR: bitset sizes do not match in copy_bits_from()" << endl;
            exit(EXIT_FAILURE);
        }

        std::memcpy(bits, other.bits, nwords * sizeof(uint64_t));
        num_1s = other.num_1s;
    }

    ~MyBitset(){
        free(bits);
    }

};

#endif