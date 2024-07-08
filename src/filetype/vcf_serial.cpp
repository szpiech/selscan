#include "vcf_serial.h"
#include "../selscan-cli.h"

#include <sstream>
#include <algorithm>
#include<cmath>
#include <set>
#include <omp.h>
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

    std::string line;
    int locus = 0;
    while (std::getline(chunk_file, line)) {
        if(line[0] == '#'){
            continue;
        }

        std::istringstream ss(line);
        std::string word;
        int word_count = 0;
        
        
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
            hapData->nhaps = (hapData->unphased)? (word_count - genotype_start_column): (word_count - genotype_start_column)*2;
        }

        if(hapData->SKIP){
            int num1s =  positions[locus].size();
            if( num1s*1.0/hapData->nhaps < hapData->MAF || 1-num1s*1.0/hapData->nhaps < hapData->MAF ){
                skiplist.push(locus);
                vector<unsigned int>().swap(positions[locus]); //clear the vector
             }  
        }
        
        locus++;
    }
    chunk_file.close();

    if(hapData->DEBUG_FLAG=="VCF") cout<<n_meta_lines<<" "<<hapData->nloci<<" "<<hapData->nhaps<<" b4 "<< locus<<" after "<<  locus-skiplist.size()<<endl; // Print results

    hapData->nloci =  (hapData->SKIP) ? locus-skiplist.size(): locus;
    hapData->initHapData(hapData->nhaps, hapData->nloci); // skip done in initHapData
    hapData->skipQueue = skiplist;
}

void VCFSerialReader::populate_positions_skipqueue() {
    //IDEA::::: MINMAF MUST BE GREATER THAN  0 that means no monomorphic site, that means all 0 sites can be considered to be filtered.
    //this way you avoid the queueu

    positions = new vector<unsigned int>[hapData->nloci];
    if(hapData->unphased)
        positions2 = new vector<unsigned int> [hapData->nloci];

    cout<<"Starting "<<"populate_positions_skipqueue ..."<<endl;
    auto start = std::chrono::high_resolution_clock::now();
    
    if(hapData->SKIP) { //prefilter all sites < MAF
        cerr << ARG_SKIP << " set. Removing all variants < " << hapData->MAF << ".\n";
        (*flog)  << ARG_SKIP << " set. Removing all variants < " << hapData->MAF << ".\n";
    }else{
        cerr << ARG_KEEP << " set. Not removing variants < " << hapData->MAF << ".\n";
        (*flog) << ARG_KEEP << " set. Not removing variants < " << hapData->MAF << ".\n";
    }

    populate_positions_skipqueue_process_chunk();

    int skipcount = 0;
    if(hapData->SKIP){
        skipcount = skiplist.size();
    }

    if(hapData->DEBUG_FLAG=="VCF") cout<<"skipcount "<<skipcount<<endl;

   if(!hapData->SKIP){
        for(int i = 0; i<= hapData->nloci-1; i++){ // i: locus_before_filter
            hapData->hapEntries[i].positions = positions[i];   
        }
    }else{
        int loc_after_skip = 0;
        for(int i = 0; i<= hapData->nloci-1; i++){ // i: locus_before_filter
            if(!skiplist.empty()){
                if(skiplist.front() == i){
                    skiplist.pop();
                    if(positions[i].size() > 0){ //skip position must have 0 length vector as it was cleared to save space
                        {
                            cout<<"ERROR: skiplist not working"<<endl;
                            throw 0;
                        }
                    }
                    continue;
                }  
            }
            hapData->hapEntries[loc_after_skip++].positions = positions[i];   
        }
    }
    
    delete[] positions;
    if(hapData->unphased){
        delete[] positions2;
    }
    //queue<int>().swap(skiplist); //clear the vector

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Serial populate_positions_skipqueue took " << duration.count() << " s." << std::endl;
}

// Function to compute symmetric difference between two vectors
void VCFSerialReader::symmetric_difference(const std::vector<unsigned int>& vec1, const std::vector<unsigned int>& vec2,  int i) {
    std::vector<unsigned int> diff;
    std::set<unsigned int> set1(vec1.begin(), vec1.end());
    std::set<unsigned int> set2(vec2.begin(), vec2.end());

    std::set_symmetric_difference(set1.begin(), set1.end(),
                                  set2.begin(), set2.end(),
                                  std::back_inserter(diff));
    hapData->hapEntries[i].xors = diff;                              
}
void VCFSerialReader::do_xor_process_chunk(int start_line, int thread_id){
    int end_line = (thread_id == num_threads - 1) ? hapData->nloci-1 : start_line + chunk_size - 1;
    for(int i = start_line; i<= hapData->nloci-1; i++){ // i: locus_after_filter
        symmetric_difference(hapData->hapEntries[i].positions, hapData->hapEntries[i-1].positions, i);
    }
}    
void VCFSerialReader::do_xor(){
    auto start = std::chrono::high_resolution_clock::now();
    hapData->hapEntries[0].xors = hapData->hapEntries[0].positions;

    // Print the number of threads
    // #pragma omp parallel
    // {
    //     #pragma omp single
    //     {
    //         std::cout << "Number of threads: " << omp_get_num_threads() << std::endl;
    //     }
    // }

    // Parallelize the computation using OpenMP
    /*
    omp_set_num_threads(this->num_threads);
    #pragma omp parallel for
    for(int i = 1; i<= hapData->nloci-1; i++){ // i: locus_after_filter
        symmetric_difference(hapData->hapEntries[i].positions, hapData->hapEntries[i-1].positions, i);
    }
    */

    threads.clear();
    
    this->chunk_size = ceil( hapData->nloci * 1.0 / num_threads); //recompute chunk size
    // if(hapData->DEBUG_FLAG=="VCF") cout<<"recomputed chunk_size "<<chunk_size<<endl; 

    for (int i = 0; i < num_threads; ++i) { // Start threads to process chunks
        int start_line = i * chunk_size; //int end_line = (i == num_threads - 1) ? num_data_lines-1 : start_line + chunk_size - 1;
        threads.emplace_back(&VCFSerialReader::do_xor_process_chunk, this, start_line, i);
    }
    for (auto& t : threads) { // Join threads
        t.join();
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Parallel xor took " << duration.count() << " s." << std::endl;
}
