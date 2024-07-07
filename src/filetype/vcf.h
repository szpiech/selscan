#ifndef VCF_PARALLEL_READER_H
#define VCF_PARALLEL_READER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
//for gzip
// #include <zlib.h>
// #include "gzstream.h"
#include<chrono>
#include "../hapmap/hapdata.h"
#include "data_reader.h"


class VCFSerialReader : public DataReader {
public:
    VCFSerialReader(const std::string filename, HapData& hapData) : DataReader(filename, hapData) {}
    /// Assigns value of HapData.nhaps and HapData.nloci  
    void get_nloci_nhaps() override {
        
    }
};

class VCFParallelReader  : public DataReader {
public:
    //VCFParallelReader(std::string& filename, HapData& hapData);
    VCFParallelReader(const std::string filename, HapData& hapData) : DataReader(filename, hapData){}
    ~VCFParallelReader(){
        //delete[] nloci_per_thread;
    }

    void getLineStartPositions(std::vector<std::streampos>& line_start_positions, int& num_meta_lines);
    void get_nloci_nhaps() override;
    void n1s_n2s_q(int* num1s_per_loci, int* num2s_per_loci, queue<int>& skiplist ) override;


private:
     std::vector<std::streampos> line_start_positions;
     std::streampos file_size; 
     int chunk_size; 

    // int lines_per_thread;

    int n_meta_lines = 0;
    int* nloci_per_thread;
    int* nlines_per_thread;
    //std::string filename;
    //int num_threads = 4;
    std::vector<std::thread> threads;
    std::mutex file_mutex;
    std::condition_variable cv;
    bool done = false;

    int nloci = 0;
    int nhaps = 0;

    inline int getGlobalLineId(int local_line_number, int thread_id){
        return local_line_number + thread_id * chunk_size;
    }
    
    // int getThreadIdFromLineNumber(int line_number){
    //     int nlines = nloci + n_meta_lines;
    //     return line_number / num_threads;

    //     if (num_threads <= 0 || nlines <= 0) {
    //         throw std::invalid_argument("numThreads and nlines must be positive integers.");
    //     }

    //     int lines_per_thread = nlines / num_threads;
    //     int remainder = nlines % num_threads;

    //     // Calculate the thread ID
    //     for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
    //         int start_line = thread_id * lines_per_thread + std::min(thread_id, remainder);
    //         int end_line = start_line + lines_per_thread + (thread_id < remainder ? 1 : 0);

    //         if (line_number >= start_line && line_number < end_line) {
    //             return thread_id;
    //         }
    //     }

    //     throw std::out_of_range("Line number is out of range.");
    // }


    /// Return the number of loci and haps in the VCF file given the thread id
    void get_nloci_nhaps_process_chunk(int start, int end, int thread_id);

    /// Return the number of 1s and 2s in the VCF file given the thread id
    void n1s_n2s_q_process_chunk(std::streampos start, std::streampos end, int thread_id, int* num1s_per_loci, int* num2s_per_loci, queue<int>& skiplist);
    
};



#endif 