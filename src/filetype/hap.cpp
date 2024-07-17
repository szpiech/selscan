#include "hap.h"

#include <sstream>
#include <algorithm>
#include "../selscan-cli.h"
#include <cmath>

using namespace std;

HapParallelReader::HapParallelReader(std::string filename, HapData* hapDataPtr) : filename(filename), hapData(*hapDataPtr){
    hapData = *hapDataPtr;
    this->filename = filename;
    flog = hapData.flog;
    num_threads = thread::hardware_concurrency(); 
  
    
    std::ifstream ifs(filename);
    ifs.seekg(0, std::ios::end);
    this->file_size = ifs.tellg();
    ifs.close();

    init_based_on_lines(); 
    populate_positions_skipqueue();
    do_xor();
}

void HapParallelReader::get_line_start_positions(std::vector<std::streampos>& start_positions){
    cout << "Running with number of concurrent threads supported: "<< num_threads << endl; 
    auto start = std::chrono::high_resolution_clock::now();
   
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file" << std::endl;
        throw 0;
    }
    std::string line;

    int num_lines = 0;    
    while (std::getline(file, line)) {
        num_lines++;
        start_positions.push_back(file.tellg() - std::streampos(line.size()) - 1); // Return end of file position if no newline found
    }
    
    file.close();

    hapData.nloci = num_lines;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;


    if(hapData.DEBUG_FLAG=="VCF") cout<<"Serial read (to get line starts and nloci) took "<< duration.count() << " s." <<endl;
}

void HapParallelReader::init_based_on_lines() {
    vector<std::streampos> all_line_start_positions;
    get_line_start_positions(all_line_start_positions);

    auto start = std::chrono::high_resolution_clock::now();
    cout<<"Starting init based on nloci "<< filename<<endl;

    int num_data_lines = all_line_start_positions.size();
    if(hapData.DEBUG_FLAG=="VCF") cout<<"num_data_lines "<<num_data_lines<<endl;

    this->chunk_size = ceil(num_data_lines * 1.0 / num_threads);
    if(hapData.DEBUG_FLAG=="VCF") cout<<"chunk_size "<<chunk_size<<endl;

    this->line_start_positions.resize(num_threads);

    for (int i = 0; i < num_threads; ++i) {
        int start_line = i * chunk_size;
        this->line_start_positions[i] = all_line_start_positions[start_line];
    }
    vector<std::streampos>().swap(all_line_start_positions); //clear the vector

    if(hapData.DEBUG_FLAG=="VCF") cout<<" "<<" Loci "<<hapData.nloci<<" "<<endl; // Print results

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Serial simple init took " << duration.count() << " s." << std::endl;    
}

void HapParallelReader::populate_positions_skipqueue_process_chunk(int start_line, int thread_id, vector<vector<vector<unsigned int> > >& positions_per_thread){
    std::ifstream chunk_file(filename);
    chunk_file.seekg(line_start_positions[thread_id]);  // go to the start of the chunk
    std::streampos endpos = (thread_id == num_threads-1) ? file_size : line_start_positions[thread_id+1];
    
    int local_line_id = 0;
    positions_per_thread[thread_id].resize(chunk_size);

    std::string line;
    while (std::getline(chunk_file, line) && chunk_file.tellg() <= endpos) {
        std::istringstream ss(line);
        std::string word;
        int word_count = 0;
        int locus = this->getGlobalLineId(local_line_id, thread_id);

        while (ss >> word) {
            ++word_count;

            int hapId = word_count - 1;
            
            char &allele1 = word[0]; 
            if (allele1 != '0' && allele1 != '1')
            {
                std::cerr << "ERROR: Alleles must be coded 0/1 only.\n";
                std::cerr << allele1 << std::endl;
                throw 0;
            }

            if(hapData.unphased){
                // char allele = '0';
                // if (allele1 == '1' && allele2 == '1'){ //allele = '2';
                //     positions2_per_thread[thread_id][local_line_id].push_back(hapId);
                // }
                // else if (allele1 == '1' || allele2 == '1'){
                //     positions_per_thread[thread_id][local_line_id].push_back(hapId);
                // }
            }
            else{
                if(allele1 == '1'){
                    positions_per_thread[thread_id][local_line_id].push_back(hapId);
                }
            }
        }
        if(local_line_id == 0){
            hapData.nhaps = hapData.unphased?(word_count/2):(word_count);
        }

        if(hapData.SKIP){
            int num1s = positions_per_thread[thread_id][local_line_id].size();
            if( num1s*1.0/hapData.nhaps < hapData.MAF || 1-num1s*1.0/hapData.nhaps < hapData.MAF ){
                skiplist_per_thread[thread_id].push(locus);
                vector<unsigned int>().swap(positions_per_thread[thread_id][local_line_id]); //clear the vector
             }  
        }
        
        local_line_id++;
    }
    chunk_file.close();
}

void HapParallelReader::populate_positions_skipqueue() {
    
    cout<<"Starting "<<"populate_positions_skipqueue ..."<<endl;
    auto start = std::chrono::high_resolution_clock::now();
    
    if(hapData.SKIP) { //prefilter all sites < MAF
        skiplist_per_thread = new queue<int>[num_threads];
        cerr << ARG_SKIP << " set. Removing all variants < " << hapData.MAF << ".\n";
        (*flog)  << ARG_SKIP << " set. Removing all variants < " << hapData.MAF << ".\n";
    }else{
        cerr << ARG_KEEP << " set. Not removing variants < " << hapData.MAF << ".\n";
        (*flog) << ARG_KEEP << " set. Not removing variants < " << hapData.MAF << ".\n";
    }
    
    positions_per_thread.resize(num_threads);
    
    if(hapData.unphased){
        positions2_per_thread = new vector<vector<unsigned int> > [num_threads];
    }

    for (int i = 0; i < num_threads; ++i) { // Start threads to process chunks
        int start_line = i * chunk_size;
        threads.emplace_back(&HapParallelReader::populate_positions_skipqueue_process_chunk, this, start_line, i, std::ref(positions_per_thread));
    }

    for (auto& t : threads) { // Join threads
        t.join();
    }

    int skipcount = 0;
    if(hapData.SKIP){
        for (int i = 0; i < num_threads; ++i) { // Aggregate results
            skipcount += skiplist_per_thread[i].size(); 
        }
    }

    if(hapData.DEBUG_FLAG=="VCF") cout<<"skipcount "<<skipcount<<endl;
    
    int num_loci_before_filter = hapData.nloci;
    hapData.nloci -= skipcount;
    int num_loci_after_filter = hapData.nloci;
    
    if(hapData.DEBUG_FLAG=="VCF") cout<<hapData.nloci<<" "<<hapData.nhaps<<" b4 "<< num_loci_before_filter<<" after "<< num_loci_after_filter<<endl; // Print results

    hapData.initHapData(hapData.nhaps, hapData.nloci); // skip done in initHapData
    hapData.skipQueue = queue<int>();


    int loc_after_skip = 0;
    for(int i = 0;  i<num_loci_before_filter; i++){ // i = loc_before_skip // do this serially
        int tid = i / chunk_size;
        int local_line_id = i % chunk_size;

        if(hapData.SKIP){
            if(!skiplist_per_thread[tid].empty()){
                if(skiplist_per_thread[tid].front() == i){
                    hapData.skipQueue.push(i);
                    skiplist_per_thread[tid].pop();
                    if(positions_per_thread[tid][local_line_id].size() > 0){
                        {
                            std::lock_guard<std::mutex> lock(file_mutex);
                            cout<<"ERROR: skiplist not working"<<endl;
                            throw 0;
                        }
                        //vector<int>().swap(positions_per_thread[tid][i]); //clear the vector
                    }
                    continue;
                }
            
            }
        }
        
        hapData.hapEntries[loc_after_skip++].positions = positions_per_thread[tid][local_line_id];            
    }

    if(hapData.SKIP){
        delete[] skiplist_per_thread;
    }
    //delete[] positions_per_thread;
    if(hapData.unphased){
        delete[] positions2_per_thread;
    }

    vector<thread>().swap(threads); //clear the vector
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Parallel populate_positions_skipqueue took " << duration.count() << " s." << std::endl;
}


void HapParallelReader::do_xor_process_chunk(int start_line, int thread_id) {
    int endline = (thread_id == num_threads-1) ? hapData.nloci-1 : start_line + chunk_size - 1;
    
    for(int locus_after_filter = start_line; locus_after_filter<=endline; locus_after_filter++){
        if(locus_after_filter==0){
            if(hapData.benchmark_flag == "XOR"){
                hapData.hapEntries[locus_after_filter].xors = hapData.hapEntries[locus_after_filter].positions;
                //hapData.hapEntries[locus_after_filter].xors1 = hapData.hapEntries[locus_after_filter].positions;
                //hapData.hapEntries[locus_after_filter].xors2 = hapData.hapEntries[locus_after_filter].positions2;
            }
        }else{
            if(hapData.benchmark_flag == "XOR"){
                vector<unsigned int>& curr_xor = hapData.hapEntries[locus_after_filter].xors;
                vector<unsigned int>& curr_positions = hapData.hapEntries[locus_after_filter].positions;
                vector<unsigned int>& prev_positions = hapData.hapEntries[locus_after_filter-1].positions;
                std::set_symmetric_difference(curr_positions.begin(), curr_positions.end(),prev_positions.begin(), prev_positions.end(),
                                std::back_inserter(curr_xor));
            }
        }
    }
}

void HapParallelReader::do_xor(){
    auto start = std::chrono::high_resolution_clock::now();

    //recompute chunk size
    this->chunk_size = ceil( hapData.nloci * 1.0 / num_threads);
    if(hapData.DEBUG_FLAG=="VCF") cout<<"recomputed chunk_size "<<chunk_size<<endl; 


    vector<thread>().swap(threads); //clear the vector
    
    for (int i = 0; i < num_threads; ++i) { // Start threads to process chunks
        int start_line = i * chunk_size; //int end_line = (i == num_threads - 1) ? num_data_lines-1 : start_line + chunk_size - 1;
        threads.emplace_back(&HapParallelReader::do_xor_process_chunk, this, start_line, i);
    }

    for (auto& t : threads) { // Join threads
        t.join();
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Parallel xor took " << duration.count() << " s." << std::endl;
}
