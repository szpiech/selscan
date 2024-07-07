#include <sstream>
using namespace std;

/****  SERIAL *****/

/****  PARALLEL *****/


/*
void VCFParallelReader::getLineStartPositions(std::vector<std::streampos>& start_positions, int& num_meta_lines){
    num_meta_lines = 0;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file" << std::endl;
        throw 0;
    }
    std::string line;
    // Read lines and record positions of newlines
    // Read the next line, which moves the file pointer to the position after the newline
    
    bool found_non_comment_line = false;
    while (std::getline(file, line)) {
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
                num_meta_lines++;
            }

        }
    }
    file.close();
}

void VCFParallelReader::get_nloci_nhaps_process_chunk(int start_line, int end_line, int thread_id) {
    nloci_per_thread[thread_id] = 0;
    std::ifstream chunk_file(filename);
     
    chunk_file.seekg(line_start_positions[thread_id]);

    std::string line;
    const int genotype_start_column = 9;
    while (std::getline(chunk_file, line) && chunk_file.tellg() <= line_start_positions[thread_id+1]-std::streampos(1)) {
        // Process the line
        nlines_per_thread[thread_id]++;

        if (line[0] == '#') {
            continue; // Skip headers
        }
    
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
        
        // Here you would parse the VCF line
        std::lock_guard<std::mutex> lock(file_mutex);
    }
    chunk_file.close();
}

void VCFParallelReader::get_nloci_nhaps() {
    done = false;

    int num_meta_lines = 0;
    vector<std::streampos> line_start_positions;
    getLineStartPositions(line_start_positions, num_meta_lines);

    auto start = std::chrono::high_resolution_clock::now();
    
    ifstream file;
    file.open(filename, std::ios::in);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }
    cout<<"Reading to get nloci nhaps "<< filename<<endl;
    nloci_per_thread = new int[num_threads];
    nlines_per_thread = new int[num_threads];


    int num_data_lines = line_start_positions.size();
    // Calculate chunk size
    int chunk_size = num_data_lines / num_threads;
    this->chunk_size = chunk_size;  
    this->line_start_positions.resize(num_threads);
    for (int i = 0; i < num_threads; ++i) {
        int start_line = i * chunk_size;
        //int end_line = (i == num_threads - 1) ? num_data_lines : start_line + chunk_size;
        this->line_start_positions[i] = line_start_positions[start_line];
    }

    int pos_iterator = 0;
    // Start threads to process chunks
    for (int i = 0; i < num_threads; ++i) {
        // std::streampos start = i * chunk_size;
        // std::streampos end = (i == num_threads - 1) ? file_size : (i + 1) * chunk_size;

        int start_line = i * chunk_size;
        int end_line = (i == num_threads - 1) ? num_data_lines : start_line + chunk_size;
        threads.emplace_back(&VCFParallelReader::get_nloci_nhaps_process_chunk, this, start_line, end_line, i);
    }

    // Join threads
    for (auto& t : threads) {
        t.join();
    }

    for (int i = 0; i < num_threads; ++i) {
        nloci += nloci_per_thread[i];
        n_meta_lines += nlines_per_thread[i];
    }
    n_meta_lines -= nloci;

    cout<<nloci<<" "<<nhaps<<endl;
    hapData.nloci = nloci;
    hapData.nhaps = nhaps;


    delete[] nloci_per_thread;
    delete[] nlines_per_thread;
    file.close();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "read nloci nhaps took " << duration.count() << " s." << std::endl;
}
*/



/*
void VCFParallelReader::n1s_n2s_q_process_chunk(std::streampos start, std::streampos end, int thread_id, int* num1s_per_loci, int* num2s_per_loci, queue<int>& skiplist){
    num1s_per_loci[hapData.nloci] = 0;
    num2s_per_loci[hapData.nloci] = 0;

    std::ifstream chunk_file(filename);
    chunk_file.seekg(start);

    std::string line;
    int num_meta_seen = 0;
    const int genotype_start_column = 9;
    int locus = 0;
    int thread_id_first_line = getThreadIdFromLineNumber(num_meta_data_lines);
    while (std::getline(chunk_file, line) && chunk_file.tellg() <= end) {
        // Process the line
        
        if (line[0] == '#') {
            num_meta_seen++;
            continue; // Skip headers
        }
        if(thread_id_first_line == thread_id){
            //prev line was meta line
            locus = 0;
        }

        // Here you would parse the VCF line
        std::lock_guard<std::mutex> lock(file_mutex);
    
        nloci_per_thread[thread_id]++;
        int startline = thread_id*nlines_per_thread[0] + 1;
        int locus = startline - num_meta_seen;
        if(thread_id != 0){
            std::istringstream ss(line);
            std::string word;
            int word_count = 0;
            while (ss >> word) {
                ++word_count;
                if(word_count<= genotype_start_column){ //TODO: check if this is correct
                    continue;
                }
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
                        num2s_per_loci[locus]++;
                    }
                    else if (allele1 == '1' || allele2 == '1'){
                        num1s_per_loci[locus]++;
                    }
                }
                else{
                    if(allele1 == '1'){
                        num1s_per_loci[locus]++;
                    }
                    if(allele2 == '1'){
                        num1s_per_loci[locus]++;
                    }
                }
            }
            //nhaps = word_count - genotype_start_column;
        }
    }
    chunk_file.close();
}
    

void VCFParallelReader::n1s_n2s_q(int* number_of_1s_per_loci, int* number_of_2s_per_loci, queue<int>& skiplist){ //can do parallely in vcf
    done = false;
    auto start = std::chrono::high_resolution_clock::now();
    
    

    ifstream file;
    file.open(filename, std::ios::in);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }
    cout<<"Reading to get n1s_n2s_q "<< filename<<endl;

    nlines_per_thread = new int[num_threads];

    //  file size,  chunk size already known
    //  Start threads to process chunks
    for (int i = 0; i < num_threads; ++i) {
        std::streampos start = i * chunk_size;
        std::streampos end = (i == num_threads - 1) ? file_size : start + chunk_size;
        threads.emplace_back(&VCFParallelReader::n1s_n2s_q_process_chunk, this, start, end, i, number_of_1s_per_loci,number_of_2s_per_loci, std::ref(skiplist));
    }

    // Join threads
    for (auto& t : threads) {
        t.join();
    }

    for (int i = 0; i < num_threads; ++i) {
        nloci += nloci_per_thread[i];
        n_meta_lines += nlines_per_thread[i];
    }
    n_meta_lines -= nloci;

    cout<<nloci<<" "<<nhaps<<endl;
    hapData.nloci = nloci;
    hapData.nhaps = nhaps;


    delete[] nloci_per_thread;
    delete[] nlines_per_thread;
    file.close();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "read nloci nhaps took " << duration.count() << " s." << std::endl;

    
    ifstream fin;
    bool unphased = hapData.unphased;

    // Pass 1: Counting so that inititalization is smooth
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());
    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    int numMapCols = 9;
    string line;
    unsigned int nloci_before_filtering = 0;
    int previous_nhaps = -1;
    int current_nhaps = 0;

    int skipcount = 0;

    int num_meta_data_lines = 0;
    while (getline(fin, line))  //Counts number of haps (cols) and number of loci (rows)
    {
        if (line[0] == '#') {
            num_meta_data_lines++;
            continue;
        }
        nloci_before_filtering++;
        current_nhaps = countFields(line) - numMapCols;


        string junk;
        char allele1, allele2, separator;
        std::stringstream ss(line);
        int number_of_1s = 0;
        int number_of_2s = 0;

        for (int i = 0; i < numMapCols; i++) {
            ss >> junk;
        }
        
        //int derived_allele_count = (unphased? (number_of_1s + number_of_2s*2) : number_of_1s);
    }
    fin.clear(); 
    fin.close();
}

*/

/*
void VCFParallelReader::n1s_n2s_q(int* num1s_per_loci, int* num2s_per_loci, queue<int>& skiplist){ //can do parallely in vcf
    bool SKIP = hapData.SKIP;
    double MAF = hapData.MAF;
    //string ARG_SKIP = hapData.ARG_SKIP; 
    string benchmark_flag = hapData.benchmark_flag;
    bool unphased = hapData.unphased;
    ofstream *flog = hapData.flog;
        
    ifstream fin;
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    int nhaps = hapData.nhaps;
    int nloci_before_filtering = hapData.nloci;
    cerr << "Loading " << nhaps << " haplotypes and " << nloci_before_filtering << " loci...\n";
    if(SKIP){
        cerr << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
        (*flog)  << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
    }

    string junk;
    char allele1, allele2, separator;
    bool skipLine = false; // to skip metadata lines

        initHapData(nhaps, nloci_before_filtering-skipcount);

        skipQueue = skiplist; 
        unsigned int nloci_after_filtering = 0;

        string prev_loc_str = "";
        string curr_loc_str = "";
        for (unsigned int locus = 0; locus < nloci_before_filtering; locus++)
        {
            curr_loc_str = "";
            for (int i = 0; i < numMapCols; i++) {
                fin >> junk;
                if (i == 0 && junk[0] == '#') { // to skip metadata lines
                    skipLine = true;
                    break;
                }
            }

            if (skipLine) { // to skip metadata lines
                getline(fin, junk);
                skipLine = false;
                locus--;
                continue;
            }

            if(!skiplist.empty()){
                if(skiplist.front() == locus){
                    skiplist.pop();    
                    getline(fin, junk);
                    continue;
                }
            }
            
            if(unphased){
                if(benchmark_flag == "XOR"){
                    hapEntries[nloci_after_filtering].xors.reserve(number_of_1s_per_loci[nloci_after_filtering]+number_of_2s_per_loci[nloci_after_filtering]);
                }
                //hapEntries[nloci_after_filtering].positions.reserve(number_of_1s_per_loci[nloci_after_filtering]+number_of_2s_per_loci[nloci_after_filtering]);

                hapEntries[nloci_after_filtering].positions.reserve(number_of_1s_per_loci[nloci_after_filtering]);
                hapEntries[nloci_after_filtering].positions2.reserve(number_of_2s_per_loci[nloci_after_filtering]);
            
            }else{
                if(benchmark_flag == "XOR"){
                    hapEntries[nloci_after_filtering].xors.reserve(number_of_1s_per_loci[nloci_after_filtering]);
                }
                hapEntries[nloci_after_filtering].positions.reserve(number_of_1s_per_loci[nloci_after_filtering]);
                

                // if(number_of_1s_per_loci[nloci_after_filtering] > nhaps/2){
                //     hapEntries[nloci_after_filtering].flipped = true;
                //     hapEntries[nloci_after_filtering].positions.reserve(nhaps - number_of_1s_per_loci[nloci_after_filtering]);

                // }else{
                //     hapEntries[nloci_after_filtering].flipped = false;
                //     hapEntries[nloci_after_filtering].positions.reserve(number_of_1s_per_loci[nloci_after_filtering]);
                // }
                
                
            }

            for (int field = 0; field <  current_nhaps ; field++)
            {
                fin >> junk;
                allele1 = junk[0];
                separator = junk[1];

                allele2 = junk[2];
                if ( (allele1 != '0' && allele1 != '1') || (allele2 != '0' && allele2 != '1') )
                {
                    cerr << "ERROR: Alleles must be coded 0/1 only.\n";
                    cerr << allele1 << " " << allele2 << endl;
                    throw 0;
                }

                //if(separator != '|'){
                //    cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
                //    throw 0;
                //}
                if(unphased){
                    char allele = '0';
                    if (allele1 == '1' && allele2 == '1'){
                        hapEntries[nloci_after_filtering].positions2.push_back(field);
                        //hapEntries[nloci_after_filtering].positions.push_back(2*field); //10
                        hapEntries[nloci_after_filtering].count2++;
                        curr_loc_str += "2";
                    }
                    else if (allele1 == '1' || allele2 == '1'){
                        hapEntries[nloci_after_filtering].positions.push_back(field);
                        //hapEntries[nloci_after_filtering].positions.push_back(2*field);
                        hapEntries[nloci_after_filtering].count1++;
                        //hapEntries[nloci_after_filtering].positions.push_back(2*field+1); //01
                        curr_loc_str += "1";

                    }else{
                        // allele = '0' implied;
                        curr_loc_str += "0";
                    }
                }
                else{
                    if(benchmark_flag != "BITSET"){
                        // if (hapEntries[nloci_after_filtering].flipped){
                        //     if(allele1 == '0'){
                        //         hapEntries[nloci_after_filtering].positions.push_back(2 * field);
                        //     }
                        //     if(allele2 == '0'){
                        //         hapEntries[nloci_after_filtering].positions.push_back(2 * field + 1);
                        //     }
                        // }else{
                        if(allele1 == '1'){
                            hapEntries[nloci_after_filtering].positions.push_back(2 * field);
                        }
                        if(allele2 == '1'){
                                hapEntries[nloci_after_filtering].positions.push_back(2 * field + 1);
                        }
                        //}
                    }
                }
            }
    }
    */