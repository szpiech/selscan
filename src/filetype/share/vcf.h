#ifndef VCF_PARALLEL_READER_H
#define VCF_PARALLEL_READER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <mutex>
#include <chrono>

using namespace std;

class VCFParallelReader {
public:
    VCFParallelReader(std::string filename, int num_threads);
    ~VCFParallelReader(){}

    /// @brief assigns value to hapData.nloci
    void get_line_start_positions(std::vector<std::streampos>& line_start_positions);
    
    void init_based_on_lines();

    // JUST AN EXAMPLE TO TEST THAT MULTITHREAING WORKS
    void get_nloci_nhaps_process_chunk(int start_line, int thread_id);
    void get_nloci_nhaps();

private:
    string filename;
    int num_threads;

    std::vector<std::thread> threads;
    std::mutex file_mutex;

    const int genotype_start_column = 9;
    int n_meta_lines = 0;

    int nloci, nhaps;

    std::vector<std::streampos> line_start_positions;
    std::streampos file_size; 
    int chunk_size; // number of lines per thread
    
    int* nloci_per_thread;  // here replace nloci with any other variable you want to calculate parallelly

    inline int getGlobalLineId(int local_line_number, int thread_id){
        return local_line_number + thread_id * this->chunk_size;
        //g = l + t * c
        //l = g % c
        //t = g / c
    }
};


#endif 