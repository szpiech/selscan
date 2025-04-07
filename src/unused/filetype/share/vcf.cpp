#include "vcf.h"

#include <sstream>
#include <algorithm>
#include<cmath>
using namespace std;
 VCFParallelReader::VCFParallelReader(std::string filename, int num_threads) {
    this->filename = filename;
    this->num_threads = num_threads;
    std::ifstream ifs(filename);
    ifs.seekg(0, std::ios::end);
    this->file_size = ifs.tellg();
    ifs.close();
    init_based_on_lines(); 
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

    nloci = num_lines - this->n_meta_lines;
    cout<<"Number of metadata lines in vcf "<< this->n_meta_lines << " lines "<< num_lines <<" nloci "<<nloci<<endl;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    cout<<"Serial read (to get line starts and nloci) took "<< duration.count() << " s." <<endl;
}

void VCFParallelReader::init_based_on_lines() {
    vector<std::streampos> all_line_start_positions;
    get_line_start_positions(all_line_start_positions);

    auto start = std::chrono::high_resolution_clock::now();
    cout<<"Starting init based on nloci "<< filename<<endl;

    int num_data_lines = all_line_start_positions.size();
    cout<<"num_data_lines "<<num_data_lines<<endl;

    this->chunk_size = ceil(num_data_lines * 1.0 / num_threads);
    cout<<"chunk_size "<<chunk_size<<endl;

    this->line_start_positions.resize(num_threads);

    for (int i = 0; i < num_threads; ++i) {
        int start_line = i * chunk_size;
        this->line_start_positions[i] = all_line_start_positions[start_line];
    }
    vector<std::streampos>().swap(all_line_start_positions); //clear the vector

    cout<<"Metalines "<<this->n_meta_lines<<" "<<" Loci "<<nloci<<" "<<endl; // Print results

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Serial simple init took " << duration.count() << " s." << std::endl;    
}


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
    while (std::getline(chunk_file, line) && chunk_file.tellg() <= endpos){ //if (line[0] == '#') continue; // Skip headers is redundant here
        nloci_per_thread[thread_id]++;
        
        if(thread_id == num_threads-1){
            std::istringstream ss(line);
            std::string word;
            int word_count = 0;
            while (ss >> word) {
                ++word_count;
            }
            nhaps = word_count - genotype_start_column;
        }
        
        
        std::lock_guard<std::mutex> lock(file_mutex);
        cout<<line.substr(6,10)<<" ";
    }
    chunk_file.close();
}

void VCFParallelReader::get_nloci_nhaps() {
    vector<std::streampos> all_line_start_positions;
    get_line_start_positions(all_line_start_positions);

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

    delete[] nloci_per_thread;  // Cleanup
    file.close();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Parallely read nloci nhaps took " << duration.count() << " s." << std::endl;

    vector<thread>().swap(threads); //clear the vector
     
}

int main(){
    VCFParallelReader vcf("test.vcf", 4);
    vcf.get_nloci_nhaps();

    return 0;
}