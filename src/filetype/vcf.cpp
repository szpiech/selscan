#include "vcf.h"

#include <sstream>
#include <algorithm>
#include "../selscan-cli.h"
#include <cmath>


using namespace std;
 VCFParallelReader::VCFParallelReader(std::string filename, HapData* hapDataPtr) : filename(filename), hapData(*hapDataPtr){
    hapData = *hapDataPtr;
    this->filename = filename;
    flog = hapData.flog;
    num_threads = hapData.num_threads;
    std::ifstream ifs(filename);
    ifs.seekg(0, std::ios::end);
    this->file_size = ifs.tellg();
    ifs.close();

    init_based_on_lines(); 
    populate_positions_skipqueue();
    do_xor();
 }

void VCFParallelReader::get_line_start_positions(std::vector<std::streampos>& start_positions){
    auto start = std::chrono::high_resolution_clock::now();
   
    this->n_meta_lines = 0;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file" << std::endl;
        throw 0;
    }
    std::string line;

    int num_lines = 0;    
    bool found_non_comment_line = false;
    while (std::getline(file, line)) {
        num_lines++;
        if(found_non_comment_line){
            start_positions.push_back(file.tellg() - std::streampos(line.size()) - 1); // Return end of file position if no newline found
        }else{
            if ( line[0] != '#') {
                if(!line.empty()){
                    std::streampos endpos = file.tellg();
                    start_positions.push_back(endpos - std::streampos(line.size() - 1)); // return the beginning of the first non-comment line
                    found_non_comment_line = true;
                }
            }else{
                this->n_meta_lines++;
            }
        }
    }
    
    file.close();

    hapData.nloci = num_lines - this->n_meta_lines;
    if(hapData.DEBUG_FLAG=="VCF") cout<<"Number of metadata lines in vcf "<< this->n_meta_lines << " lines "<< num_lines <<" nloci "<<hapData.nloci<<endl;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    if(hapData.DEBUG_FLAG=="VCF") cout<<"Serial read (to get line starts and nloci) took "<< duration.count() << " s." <<endl;
}

void VCFParallelReader::init_based_on_lines() {
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

    if(hapData.DEBUG_FLAG=="VCF") cout<<"Metalines "<<this->n_meta_lines<<" "<<" Loci "<<hapData.nloci<<" "<<endl; // Print results

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Serial simple init took " << duration.count() << " s." << std::endl;    
}

/*
void VCFParallelReader::get_nloci_nhaps_process_chunk(int start_line, int thread_id) {
    nloci_per_thread[thread_id] = 0;    // initialize the number of loci for this thread

    std::ifstream chunk_file(filename);
    chunk_file.seekg(line_start_positions[thread_id]);  // go to the start of the chunk

    const int genotype_start_column = 9;
    std::streampos endpos;
    if(thread_id == num_threads-1){
        endpos = file_size ; 
    }else{
        endpos = line_start_positions[thread_id+1] ;
    }
    
    std::string line;
    while (std::getline(chunk_file, line) && chunk_file.tellg() <= endpos) { //if (line[0] == '#') continue; // Skip headers is redundant here
        nloci_per_thread[thread_id]++;
        // 
        // if(thread_id == num_threads-1){
        //     std::istringstream ss(line);
        //     std::string word;
        //     int word_count = 0;
        //     while (ss >> word) {
        //         ++word_count;
        //     }
        //     nhaps = word_count - genotype_start_column;
        // }
        // 
        
        // std::lock_guard<std::mutex> lock(file_mutex);
        // cout<<line.substr(6,10)<<" ";
    }
    chunk_file.close();
}

void VCFParallelReader::get_nloci_nhaps() {
    vector<std::streampos> all_line_start_positions;
    get_line_start_positions(all_line_start_positions, this->n_meta_lines);

    auto start = std::chrono::high_resolution_clock::now();
    
    ifstream file;
    file.open(filename, std::ios::in);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }
    cout<<"Reading to get nloci nhaps "<< filename<<endl;
    nloci_per_thread = new int[num_threads];

    int num_data_lines = all_line_start_positions.size();
    cout<<"num_data_lines "<<num_data_lines<<endl;

    int chunk_size = ceil(num_data_lines * 1.0 / num_threads);
    cout<<"chunk_size "<<chunk_size<<endl;

    this->chunk_size = chunk_size;  
    this->line_start_positions.resize(num_threads);

    for (int i = 0; i < num_threads; ++i) {
        int start_line = i * chunk_size;
        this->line_start_positions[i] = all_line_start_positions[start_line];
    }

    for (int i = 0; i < num_threads; ++i) { // Start threads to process chunks
        int start_line = i * chunk_size;
        //int end_line = (i == num_threads - 1) ? num_data_lines-1 : start_line + chunk_size - 1;
        threads.emplace_back(&VCFParallelReader::get_nloci_nhaps_process_chunk, this, start_line, i);
    }

    for (auto& t : threads) { // Join threads
        t.join();
    }

    for (int i = 0; i < num_threads; ++i) { // Aggregate results
        nloci += nloci_per_thread[i];
    }

    cout<<n_meta_lines<<" "<<nloci<<" "<<nhaps<<endl; // Print results
    hapData.nloci = nloci;
    hapData.nhaps = nhaps;

    delete[] nloci_per_thread;  // Cleanup
    file.close();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Parallely read nloci nhaps took " << duration.count() << " s." << std::endl;

    vector<thread>().swap(threads); //clear the vector
     
}
*/

void VCFParallelReader::populate_positions_skipqueue_process_chunk(int start_line, int thread_id){
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
            if(word_count<= genotype_start_column){ //TODO: check if this is correct
                continue;
            }

            int hapId = word_count - genotype_start_column - 1;
            
            char &allele1 = word[0]; //separator = word[1];
            char &allele2 = word[2];
            if ( (allele1 != '0' && allele1 != '1') || (allele2 != '0' && allele2 != '1') )
            {
                std::cerr << "ERROR: Alleles must be coded 0/1 only.\n";
                std::cerr << allele1 << " " << allele2 << std::endl;
                throw 0;
            }

            if(hapData.unphased){
                char allele = '0';
                if (allele1 == '1' && allele2 == '1'){ //allele = '2';
                    positions2_per_thread[thread_id][local_line_id].push_back(hapId);
                }
                else if (allele1 == '1' || allele2 == '1'){
                    positions_per_thread[thread_id][local_line_id].push_back(hapId);
                }
            }
            else{
                if(allele1 == '1'){
                    positions_per_thread[thread_id][local_line_id].push_back(2*hapId);
                }
                if(allele2 == '1'){
                    positions_per_thread[thread_id][local_line_id].push_back(2*hapId+1);
                }
            }
        }
        if(local_line_id == 0){
            hapData.nhaps = hapData.unphased?(word_count - genotype_start_column):(word_count - genotype_start_column)*2;
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

void VCFParallelReader::populate_positions_skipqueue() {
    
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
    
    positions_per_thread = new vector<vector<unsigned int> > [num_threads];
    
    if(hapData.unphased){
        positions2_per_thread = new vector<vector<unsigned int> > [num_threads];
    }

    for (int i = 0; i < num_threads; ++i) { // Start threads to process chunks
        int start_line = i * chunk_size;
        threads.emplace_back(&VCFParallelReader::populate_positions_skipqueue_process_chunk, this, start_line, i);
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
    
    if(hapData.DEBUG_FLAG=="VCF") cout<<n_meta_lines<<" "<<hapData.nloci<<" "<<hapData.nhaps<<" b4 "<< num_loci_before_filter<<" after "<< num_loci_after_filter<<endl; // Print results

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
    delete[] positions_per_thread;
    if(hapData.unphased){
        delete[] positions2_per_thread;
    }

    vector<thread>().swap(threads); //clear the vector
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Parallel populate_positions_skipqueue took " << duration.count() << " s." << std::endl;
}


void VCFParallelReader::do_xor_process_chunk(int start_line, int thread_id) {
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

void VCFParallelReader::do_xor(){
    auto start = std::chrono::high_resolution_clock::now();

    //recompute chunk size
    this->chunk_size = ceil( hapData.nloci * 1.0 / num_threads);
    if(hapData.DEBUG_FLAG=="VCF") cout<<"recomputed chunk_size "<<chunk_size<<endl; 


    vector<thread>().swap(threads); //clear the vector
    
    for (int i = 0; i < num_threads; ++i) { // Start threads to process chunks
        int start_line = i * chunk_size; //int end_line = (i == num_threads - 1) ? num_data_lines-1 : start_line + chunk_size - 1;
        threads.emplace_back(&VCFParallelReader::do_xor_process_chunk, this, start_line, i);
    }

    for (auto& t : threads) { // Join threads
        t.join();
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Parallel xor took " << duration.count() << " s." << std::endl;
    
    //vector<thread>().swap(threads); //clear the vector
    
    /*
    for (int locus_after_filter = 0; locus_after_filter < hapData.nloci; locus_after_filter++){
        if(locus_after_filter==0){
            if(hapData.benchmark_flag == "XOR"){
                vector<unsigned int>& source = hapData.hapEntries[locus_after_filter].positions;
                vector<unsigned int>& destination = hapData.hapEntries[locus_after_filter].xors;
                std::copy(source.begin(), source.end(), destination.begin());
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
    */

    /* flipped 
    for (int locus_after_filter = 0; locus_after_filter < hapData.nloci; locus_after_filter++){
        // vector<unsigned int>& one_positions = hapData.hapEntries[locus_after_filter].positions;
        if(hapData.hapEntries[locus_after_filter].positions.size() > hapData.nhaps/2){
            hapData.hapEntries[locus_after_filter].flipped = true;

            vector<unsigned int> copy_pos;
            copy_pos.reserve(this->nhaps - hapData.hapEntries[locus_after_filter].positions.size());
            int cnt = 0;
            for(int i = 0; i< this->nhaps; i++){
                unsigned int curr =  hapData.hapEntries[locus_after_filter].positions[cnt];
                if(i==curr){
                    cnt++;
                }else{
                    copy_pos.push_back(i);
                }
            }
            
            hapData.hapEntries[locus_after_filter].positions = copy_pos;
            copy_pos.clear();
        }else{
            hapData.hapEntries[locus_after_filter].flipped = false;
        }
    }
    */

}
