#include "vcf_serial.h"
#include "../selscan-cli.h"

#include <sstream>
#include <algorithm>
#include<cmath>
using namespace std;


void VCFSerialReader::get_line_start_positions(){
    auto start = std::chrono::high_resolution_clock::now();
   
    this->n_meta_lines = 0;
    igzstream file(filename.c_str());
    if (file.fail()) {
        std::cerr << "Error opening file" << std::endl;
        throw 0;
    }

    std::string line;
    int num_lines = 0;    
    bool found_non_comment_line = false;
    while (std::getline(file, line)) {
        num_lines++;
        if(!found_non_comment_line){
            if ( line[0] != '#') {
                if(!line.empty()){
                    found_non_comment_line = true;
                }
            }else{
                this->n_meta_lines++;
            }
        }
    }
    
    file.close();

    hapData->nloci = num_lines - this->n_meta_lines;
    if(hapData->DEBUG_FLAG=="VCF") cout<<"Number of metadata lines in vcf "<< this->n_meta_lines << " lines "<< num_lines <<" nloci "<<hapData->nloci<<endl;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    if(hapData->DEBUG_FLAG=="VCF") cout<<"Serial read (to get line starts and nloci) took "<< duration.count() << " s." <<endl;
}

void VCFSerialReader::init_based_on_lines() {
    get_line_start_positions();

    auto start = std::chrono::high_resolution_clock::now();
    cout<<"Starting init based on nloci "<< filename<<endl;

    this->chunk_size = ceil(hapData->nloci * 1.0 / num_threads);
    if(hapData->DEBUG_FLAG=="VCF") cout<<"chunk_size "<<chunk_size<<endl;

    if(hapData->DEBUG_FLAG=="VCF") cout<<"Metalines "<<this->n_meta_lines<<" "<<" Loci "<<hapData->nloci<<" "<<endl; // Print results

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Serial simple init took " << duration.count() << " s." << std::endl;    
}

void VCFSerialReader::populate_positions_skipqueue_process_chunk(){
    igzstream chunk_file(filename.c_str());
    if(chunk_file.fail()){
        std::cerr << "Error opening file" << std::endl;
        throw 0;
    }

    positions = new vector<unsigned int> [hapData->nloci];


    if(hapData->unphased)
        positions2 = new vector<unsigned int> [hapData->nloci];


    std::string line;
    int locus = 0;
    while (std::getline(chunk_file, line)) {
        if(line[0] == '#'){
            continue;
        }

        std::istringstream ss(line);
        std::string word;
        int word_count = 0;
        
        positions[locus].reserve(hapData->nhaps);
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

            if(hapData->unphased){
                char allele = '0';
                if (allele1 == '1' && allele2 == '1'){ //allele = '2';
                    positions2[locus].push_back(hapId);
                }
                else if (allele1 == '1' || allele2 == '1'){
                    positions[locus].push_back(hapId);
                }
            }
            else{
                if(allele1 == '1'){
                    positions[locus].push_back(2*hapId);
                }
                if(allele2 == '1'){
                    positions[locus].push_back(2*hapId+1);
                }
            }
        }
        if(locus == 0){
            hapData->nhaps = word_count - genotype_start_column;
        }

        if(hapData->SKIP){
            int num1s = positions[locus].size();
            if( num1s*1.0/hapData->nhaps < hapData->MAF || 1-num1s*1.0/hapData->nhaps < hapData->MAF ){
                skiplist.push(locus);
                vector<unsigned int>().swap(positions[locus]); //clear the vector
             }  
        }
        
        locus++;
    }
    chunk_file.close();
}

void VCFSerialReader::populate_positions_skipqueue() {
    
    cout<<"Starting "<<"populate_positions_skipqueue ..."<<endl;
    auto start = std::chrono::high_resolution_clock::now();
    
    if(hapData->SKIP) { //prefilter all sites < MAF
        cerr << ARG_SKIP << " set. Removing all variants < " << hapData->MAF << ".\n";
        (*flog)  << ARG_SKIP << " set. Removing all variants < " << hapData->MAF << ".\n";
    }else{
        cerr << ARG_KEEP << " set. Not removing variants < " << hapData->MAF << ".\n";
        (*flog) << ARG_KEEP << " set. Not removing variants < " << hapData->MAF << ".\n";
    }

    // do serial procsessing
    populate_positions_skipqueue_process_chunk();

    int skipcount = 0;
    if(hapData->SKIP){
        skipcount = skiplist.size();
    }

    if(hapData->DEBUG_FLAG=="VCF") cout<<"skipcount "<<skipcount<<endl;
    
    int num_loci_before_filter = hapData->nloci;
    hapData->nloci -= skipcount;
    int num_loci_after_filter = hapData->nloci;
    
    if(hapData->DEBUG_FLAG=="VCF") cout<<n_meta_lines<<" "<<hapData->nloci<<" "<<hapData->nhaps<<" b4 "<< num_loci_before_filter<<" after "<< num_loci_after_filter<<endl; // Print results

    hapData->initHapData(hapData->nhaps, hapData->nloci); // skip done in initHapData
    hapData->skipQueue = skiplist;
    
    int loc_after_skip = 0;
    for(int i = 0;  i<num_loci_before_filter; i++){ // i = loc_before_skip // do this serially
        if(hapData->SKIP){
            if(!skiplist.empty()){
                if(skiplist.front() == i){
                    skiplist.pop();
                    if(positions[i].size() > 0){ //skip position must have 0 length vector as it was cleared to save space
                        {
                            cout<<"ERROR: skiplist not working"<<endl;
                            throw 0;
                        }
                        //vector<int>().swap(positions_per_thread[tid][i]); //clear the vector
                    }
                    continue;
                }  
            }
        }
        
       // hapData.hapEntries[loc_after_skip++].positions.resize(positions[i].size());
       // std::copy(positions[i].begin(), positions[i].end(), hapData.hapEntries[loc_after_skip].positions.begin());
       positions[i].shrink_to_fit();
       //hapData.hapEntries[loc_after_skip++].positions.swap(positions[i]); 

       hapData->hapEntries[loc_after_skip++].positions = (positions[i]);   

    }

    //queue<int>().swap(skiplist); //clear the vector
    delete[] positions;

    if(hapData->unphased){
        delete[] positions2;
    }
    //queue<int>().swap(skiplist); //clear the vector

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Serial populate_positions_skipqueue took " << duration.count() << " s." << std::endl;
}


void VCFSerialReader::do_xor_process_chunk(int start_line, int thread_id) {
    int endline = (thread_id == num_threads-1) ? hapData->nloci-1 : start_line + chunk_size - 1;
    cout<<"stend "<< start_line<<" "<<endline<<endl;
    for(int locus_after_filter = start_line; locus_after_filter<=endline; locus_after_filter++){

        if(locus_after_filter==0){
            if(hapData->benchmark_flag == "XOR"){
                hapData->hapEntries[locus_after_filter].xors = hapData->hapEntries[locus_after_filter].positions;
                //cout<<"p0 size"<<hapData.hapEntries[locus_after_filter].xors.size()<<endl;
                //hapData.hapEntries[locus_after_filter].xors1 = hapData.hapEntries[locus_after_filter].positions;
                //hapData.hapEntries[locus_after_filter].xors2 = hapData.hapEntries[locus_after_filter].positions2;
            }
        }else{
            if(hapData->benchmark_flag == "XOR"){
                vector<unsigned int> curr_xor;
                std::set_symmetric_difference(hapData->hapEntries[locus_after_filter].positions.begin(), hapData->hapEntries[locus_after_filter].positions.end(),hapData->hapEntries[locus_after_filter-1].positions.begin(), hapData->hapEntries[locus_after_filter-1].positions.end(),
                                std::back_inserter(curr_xor));
                hapData->hapEntries[locus_after_filter].xors.resize(curr_xor.size());
                hapData->hapEntries[locus_after_filter].xors = (curr_xor);
            }
        }
    }
}

void VCFSerialReader::do_xor(){
    auto start = std::chrono::high_resolution_clock::now();
    threads.clear();
    //recompute chunk size
    this->chunk_size = ceil( hapData->nloci * 1.0 / num_threads);
    if(hapData->DEBUG_FLAG=="VCF") cout<<"recomputed chunk_size "<<chunk_size<<endl; 

    cout<<num_threads<<" threads "<<endl;
    // for (int i = 0; i < num_threads; ++i) { // Start threads to process chunks
    //     int start_line = i * chunk_size; //int end_line = (i == num_threads - 1) ? num_data_lines-1 : start_line + chunk_size - 1;
    //     threads.emplace_back(&VCFSerialReader::do_xor_process_chunk, this, start_line, i);
    // }
    // for (auto& t : threads) { // Join threads
    //     t.join();
    // }


    // int endline = (thread_id == num_threads-1) ? hapData.nloci-1 : start_line + chunk_size - 1;
    // cout<<"stend "<< start_line<<" "<<endline<<endl;
    for(int locus_after_filter = 0; locus_after_filter<=hapData->nloci-1; locus_after_filter++){
        cout<<locus_after_filter<<endl;
        if(locus_after_filter==0){
            if(hapData->benchmark_flag == "XOR"){
                hapData->hapEntries[locus_after_filter].xors = hapData->hapEntries[locus_after_filter].positions;
                cout<<"p0 size"<<hapData->hapEntries[locus_after_filter].xors.size()<<endl;
                //hapData.hapEntries[locus_after_filter].xors1 = hapData.hapEntries[locus_after_filter].positions;
                //hapData.hapEntries[locus_after_filter].xors2 = hapData.hapEntries[locus_after_filter].positions2;
            }
        }else{
            if(hapData->benchmark_flag == "XOR"){
                        vector<unsigned int>& curr_xor = hapData->hapEntries[locus_after_filter].xors;
        vector<unsigned int>& curr_positions = hapData->hapEntries[locus_after_filter].positions;
        hapData->hapEntries[locus_after_filter].positions.shrink_to_fit();
        vector<unsigned int>& prev_positions = hapData->hapEntries[locus_after_filter-1].positions;
        hapData->hapEntries[locus_after_filter-1].positions.shrink_to_fit();

        std::set_symmetric_difference(curr_positions.begin(), curr_positions.end(),prev_positions.begin(), prev_positions.end(), std::back_inserter(curr_xor));  
                
                /*
                //vector<unsigned int> curr_xor;
                std::set_symmetric_difference(hapData.hapEntries[locus_after_filter].positions.begin(), hapData.hapEntries[locus_after_filter].positions.end(),hapData.hapEntries[locus_after_filter-1].positions.begin(), hapData.hapEntries[locus_after_filter-1].positions.end(),
                                std::back_inserter(hapData.hapEntries[locus_after_filter].xors));
                //hapData.hapEntries[locus_after_filter].xors.resize(curr_xor.size());
                //hapData.hapEntries[locus_after_filter].xors = (curr_xor);
                */
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Parallel xor took " << duration.count() << " s." << std::endl;
    vector<thread>().swap(threads); //clear the vector
    
}
