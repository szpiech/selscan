#ifndef VCF_PARALLEL_READER_H
#define VCF_PARALLEL_READER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <mutex>
#include <chrono>
#include "../hapmap/hapdata.h"

class VCFParallelReader {
public:
    VCFParallelReader(std::string filename, HapData* hapDataPtr);
    ~VCFParallelReader(){}

    /// @brief assigns value to hapData.nloci
    void get_line_start_positions(std::vector<std::streampos>& line_start_positions);
    
    void init_based_on_lines();

    /// @brief Assigns value of #hapData.positions  #hapData.positions2  #skiplist  #hapData.nhaps
    void populate_positions_skipqueue();
    void do_xor();

private:
    string filename;
    HapData& hapData;
    ofstream* flog;
    int num_threads;

    std::vector<std::thread> threads;
    std::mutex file_mutex;

    const int genotype_start_column = 9;
    int n_meta_lines = 0;

    std::vector<std::streampos> line_start_positions;
    std::streampos file_size; 
    int chunk_size; // number of lines per thread
    
    queue<int>* skiplist_per_thread;
    vector<vector<unsigned int> > * positions_per_thread;
    vector<vector<unsigned int> > * positions2_per_thread;

    inline int getGlobalLineId(int local_line_number, int thread_id){
        return local_line_number + thread_id * this->chunk_size;
        //g = l + t * c
        //l = g % c
        //t = g / c
    }

    /// Return the number of 1s and 2s in the VCF file given the thread id
    void populate_positions_skipqueue_process_chunk(int start_line, int thread_id);  
    void do_xor_process_chunk(int start_line, int thread_id);
};


#endif 