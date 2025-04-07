#ifndef __VCF_SERIAL_READER_H__
#define __VCF_SERIAL_READER_H__

//for gzip
#include <zlib.h>
#include "../gzstream.h"
#include "../hapmap/hapdata.h"
#include <thread>
#include <mutex>
#include <chrono>


class VCFSerialReader   {
public:
    std::string filename;
    HapData* hapData;
    VCFSerialReader(std::string filename, HapData* hapData);
    ~VCFSerialReader(){
        threads.clear();
        //delete[] nloci_per_thread;
    }

    /// @brief assigns value to hapData.nloci
    /// @param line_start_positions 
    void get_line_start_positions();
    
    void init_based_on_lines() ;

    /// @brief Assigns value of #hapData.positions  #hapData.positions2  #skiplist  #hapData.nhaps
    void populate_positions_skipqueue() ;
    void do_xor() ;

private:
    const int genotype_start_column = 9;
    int n_meta_lines = 0;

    std::vector<std::thread> threads;
    std::mutex file_mutex;
    int chunk_size;
    int num_threads;

    ofstream* flog;

    queue<int> skiplist;
    vector<vector<unsigned int> > positions;
    vector<unsigned int> * positions2;

    /// Return the number of 1s and 2s in the VCF file given the thread id
    void populate_positions_skipqueue_process_chunk();  
    //void do_xor_process_chunk(int start_line, int thread_id);
    //void do_xor_process_chunk();
    void do_xor_process_chunk(int start_line, int thread_id);

    void symmetric_difference(const std::vector<unsigned int>& vec1, const std::vector<unsigned int>& vec2, int i);


};


#endif