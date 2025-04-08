
#include "selscan-data.h"


#include <set>



void HapMap::initParamsInHap(HapData &hapData){
    hapData.MISSING_ALLOWED = p.MISSING_ALLOWED;
    hapData.unphased = p.UNPHASED;
    hapData.MULTI_MAF = p.MULTI_MAF;

    //hapData.MULTI_CHR = p.MULTI_CHR;
    // if(p.MULTI_CHR){
    //     cerr<<"Multi-chromosome support enabled."<<endl;
    //     // MULTI-chromosome support added
        
    //     if(p.CHR_LIST!=""){
    //         std::stringstream ss(p.CHR_LIST);
    //         std::string item;
    //         while (std::getline(ss, item, ',')) { // Split the input string by commas
    //             if(item!="")
    //                 chr_set.insert(item);
    //         }
    //     }

    //     if(chr_set.empty()){
    //         cerr<<"ERROR: No valid chromosome list provided.\n";
    //     }
    // }

    // if we are in multi param mode, or using xp, we want to set the MAF cutoff to the minimum
    this->MIN_MAF_CUTOFF = p.MAF;
    if(p.MULTI_MAF){
        for (const auto& param : ps) {
            if(param.MAF < MIN_MAF_CUTOFF){
                this->MIN_MAF_CUTOFF = param.MAF;
            }
        }
    }
}




// ms hap format : number of rows = number of haplotypes, number of columns = number of loci
void HapMap::readHapDataMSHap(string filename, HapData &hapData)
{
    //HANDLE_ERROR("Bug found in ms-style hap format. This function is undergoing development.");
    initParamsInHap(hapData);



    // ifstream infile(filename); // Replace with your actual file path
    // if (!infile) {
    //     cerr << "Failed to open file." << endl;
    //     return;
    // }

    // vector<vector<int>> raw_matrix;
    // string line;

    // // Step 1: Read file into a 2D vector
    // while (getline(infile, line)) {
    //     vector<int> row;
    //     for (char ch : line) {
    //         if (ch == '0' || ch == '1') {
    //             row.push_back(ch - '0');
    //         }
    //     }
    //     raw_matrix.push_back(row);
    // }

    // infile.close();

    // if (raw_matrix.empty()) {
    //     cerr << "Matrix is empty or invalid." << endl;
    //     return;
    // }

    // size_t num_rows = raw_matrix.size();
    // size_t num_cols = raw_matrix[0].size();

    // // Step 2: Count zeros in each column
    // vector<int> zero_counts(num_cols, 0);
    // for (size_t j = 0; j < num_cols; ++j) {
    //     for (size_t i = 0; i < num_rows; ++i) {
    //         if (raw_matrix[i][j] == 0) {
    //             zero_counts[j]++;
    //         }
    //     }
    // }

    // // Step 3: Collect indices of columns with > 5 zeros
    // vector<size_t> selected_cols;
    // for (size_t j = 0; j < num_cols; ++j) {
    //     if (zero_counts[j]/50.0 >= p.MAF && 1-zero_counts[j]/50.0 >= p.MAF) {
    //         selected_cols.push_back(j);
    //     }
    // }

    // // Step 4: Build filtered matrix
    // vector<vector<int>> filtered_matrix;
    // for (const auto& row : raw_matrix) {
    //     vector<int> new_row;
    //     for (size_t j : selected_cols) {
    //         new_row.push_back(row[j]);
    //     }
    //     filtered_matrix.push_back(new_row);
    // }

    // // Output the filtered matrix
    // cout << "Filtered Matrix (columns with >5 zeros):" << endl;
    // for (const auto& row : filtered_matrix) {
    //     cout<< row.size() << " " <<8324-row.size() << endl;
    //     break;
    //     // for (int val : row) {
    //     //     cout << val;

    //     // }
    //     // cout << endl;
    // }
    // cout << "Number of rows: " << filtered_matrix.size() << endl;

    //return ;


    // Phase === PHASE 1: Determine nloci, nhaps, and MAF-based skiplist ===
    //PHASE 1: Read input file to get "nloci", "nhaps" and "skiplist"
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        exit(EXIT_FAILURE);
    }

    string line;
    string prev_line;
    int current_nhaps = 0;

    vector<int> skiplist;

    getline(fin, prev_line); // first line
    prev_line.erase(std::remove_if(prev_line.begin(), prev_line.end(), [](char c) {
            return std::isspace(static_cast<unsigned char>(c));
        }), prev_line.end());  // erase spaces inside and outside

    int nloci_before_filter =  prev_line.length(); // number of loci before MAF filtering

    vector<int> number_of_1s_per_locus(nloci_before_filter, 0);
    vector<int> number_of_2s_per_locus(nloci_before_filter, 0);

 //Counts number of haps (rows) and number of loci (cols)
 //if any lines differ, send an error message and throw an exception

 int hap = 1;
 while (getline(fin, line))
 {
     line.erase(std::remove_if(line.begin(), line.end(), [](char c) {
         return std::isspace(static_cast<unsigned char>(c));
     }), line.end());  // erase spaces inside and outside

     if(prev_line.length() != line.length()){
         HANDLE_ERROR("ERROR: line " + std::to_string(hap) + " of " + filename + " has " + std::to_string(line.length()) +
                       " loci, but the previous line has " + std::to_string(prev_line.length()));
     }

     for(int loc=0; loc<nloci_before_filter; loc++){
         if(p.UNPHASED){
             if (hap % 2 == 1){
                 if(line[loc]=='1'  && prev_line[loc]=='1'){
                     number_of_2s_per_locus[loc]++;
                 }else if ( (line[loc]=='1' && prev_line[loc]=='0') || (line[loc]=='0' && prev_line[loc]=='1') ){ // ==1
                     number_of_1s_per_locus[loc]++;
                 }
             }
         }else{
             if (hap % 2 == 1){
                 if(line[loc]=='1'){
                     number_of_1s_per_locus[loc]++;
                 }
                 if(prev_line[loc]=='1'){
                     number_of_1s_per_locus[loc]++;
                 }
             }
             
         }
     }
         
     prev_line = line; // store the previous line
     hap++;
 }
 
 int num_rows = hap;
 current_nhaps = (p.UNPHASED) ? num_rows/2 : num_rows;

 for(int loc=0; loc<nloci_before_filter; loc++){
     if(shouldSkipLocus(number_of_1s_per_locus[loc], number_of_2s_per_locus[loc], num_rows)){
         skiplist.push_back(loc);
     }
 }
 
 fin.clear(); // clear error flags
 fin.close();

 number_of_1s_per_locus.clear();
 number_of_2s_per_locus.clear();

// === Logging summary ===

if (p.SKIP) {
    LOG(ARG_SKIP << " set. Removing all variants < " << p.MAF);
    LOG("Removed " << skiplist.size() << " low frequency variants from haplotype data.");
} else {
    LOG(ARG_KEEP << " set. NOT removing variants < " << p.MAF);
}
 
    hapData.initHapData(current_nhaps, nloci_before_filter - skiplist.size());
    
    //PHASE 2: Open input File 2nd time To Load into Data Structure
    LOG("Re-loading input to read " << current_nhaps << " haplotypes and "
        << (nloci_before_filter - skiplist.size()) << " loci into data structures.");
    
    fin.open(filename.c_str());
    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        exit(EXIT_FAILURE);
    }

    stringstream ss;
    char allele1;
    int locus_after_filter = 0;

    // getline(fin, prev_line);
    // ss.str(prev_line);
    // prev_line.erase(std::remove_if(prev_line.begin(), prev_line.end(), [](char c) {
    //         return std::isspace(static_cast<unsigned char>(c));
    //     }), prev_line.end());
   
    prev_line = "";
    for (int hap = 0; hap < num_rows; hap++)
    {
        getline(fin, line);
        line.erase(std::remove_if(line.begin(), line.end(), [](char c) {
            return std::isspace(static_cast<unsigned char>(c));
        }), line.end());   

        int skip_idx = 0;

        int loc_after_filter = 0;
        for(int loc=0; loc<nloci_before_filter; loc++){
            
            if (skip_idx < skiplist.size()) {
                if ( loc == skiplist[skip_idx]) {
                    if(hap==num_rows-1){
                        hapData.skipQueue.push(loc);
                    }
                    ++skip_idx;  // move to the next skip index
                    continue;    // skip printing this locus
                }
            }
            
            if(p.UNPHASED){
                if (hap % 2 == 1){
                    if(line[loc]=='1'  && prev_line[loc]=='1'){
                        hapData.addAllele2(loc_after_filter, (hap-1)/2);
                    }else if ( (line[loc]=='1' && prev_line[loc]=='0') || (line[loc]=='0' && prev_line[loc]=='1') ){ // ==1
                        hapData.addAllele1(loc_after_filter, (hap-1)/2);
                    }
                }
            }else{
                if (hap % 2 == 1){
                    if((prev_line[loc]!='1' && prev_line[loc]!='0') || (line[loc]!='1' && line[loc]!='0')){
                        HANDLE_ERROR("Alleles must be coded 0/1 only.");
                    }

                    if(prev_line[loc]=='1'){
                        hapData.addAllele1(loc_after_filter, (hap-1));
                    } 
                    if (line[loc]=='1'){
                        hapData.addAllele1(loc_after_filter, hap);
                    }

                }
                
            }
            loc_after_filter++;
        }

        prev_line = line; // store the previous line
        
    }
    fin.close();

}




//reads in haplotype data and also does basic checks on integrity of format
//returns a populated HaplotypeData structure if successful
//impute hap IMPUTE HAP is transposed format where row represents loci,  column replesent individual
//so wc -l of impute hap is same as wc -l of map.
void HapMap::readHapDataTHAP(string filename, HapData& hapData)
{
    initParamsInHap(hapData);

    //verification
    // ifstream infile(filename); // Input file path
    // if (!infile) {
    //     cerr << "Failed to open file." << endl;
    //     return;
    // }

    // vector<vector<int>> matrix;
    // string line;

    // int nrowss = 0;
    // int ncolss = 0;
    // // Step 1: Read the matrix from file (text format)
    // while (getline(infile, line)) {
    //     vector<int> row;
    //     ncolss = 0;
    //     for (char ch : line) {
    //         if (ch == '0' || ch == '1') {
    //             row.push_back(ch - '0');
    //             ncolss++;
    //         }
    //     }
    //     if (!row.empty()) {
    //         matrix.push_back(row);
    //     }
    //     nrowss++;
    // }

    // infile.close();

    // if (matrix.empty()) {
    //     cerr << "Matrix is empty or invalid." << endl;
    //     return;
    // }

    // // Step 2: Filter rows with more than 5 zeros
    // vector<vector<int>> filtered_matrix;
    // for (const auto& row : matrix) {
    //     int zero_count = 0;
    //     for (int val : row) {
    //         if (val == 0) zero_count++;
    //     }
    //     if (zero_count*1.0/ncolss > p.MAF && 1-zero_count*1.0/ncolss > p.MAF) {
    //         filtered_matrix.push_back(row);
    //     }
    // }

    // // Step 3: Print the filtered matrix
    // cout << "Filtered Matrix (rows with >5 zeros):" << endl;
    // for (const auto& row : filtered_matrix) {
    //     // for (int val : row) {
    //     //     cout << val;
    //     // }
    //     // cout << endl;
    // }
    // cout << "Number of rows (loci): " << filtered_matrix.size() << endl;


    //PHASE 1: Read input file to get "nloci", "nhaps" and "skiplist"
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        HANDLE_ERROR("Failed to open " + filename + " for reading.");
    }

    string line;
    
    int previous_nhaps = -1;
    int current_nhaps = 0;

    queue<int> skiplist;
    int nloci_before_filter = 0;

    //Counts number of haps (rows) and number of loci (cols)
    //if any lines differ, send an error message and throw an exception
    while (getline(fin, line))
    {
        nloci_before_filter++;
        line.erase(std::remove_if(line.begin(), line.end(), [](char c) {
            return std::isspace(static_cast<unsigned char>(c));
        }), line.end());
        current_nhaps = line.length(); // number of haps 
        
        int number_of_2s = 0; // for unphased
        int number_of_1s = 0;
        for (int hap = 0; hap < current_nhaps; hap++){ 
            if(p.UNPHASED){
                if (hap % 2 == 0){
                    if(line[hap]=='1'  && line[hap+1]=='1'){
                        number_of_2s++;
                    }else if ( (line[hap]=='1' && line[hap+1]=='0') || (line[hap]=='0' && line[hap+1]=='1') ){ // ==1
                        number_of_1s++;
                    }
                }
            }else{
                if(line[hap]=='1'){
                    number_of_1s++;
                }
            }
        }

        if(shouldSkipLocus(number_of_1s, number_of_2s, current_nhaps)){
            skiplist.push(nloci_before_filter-1);
        }
        
    
        if (previous_nhaps < 0)
        {
            previous_nhaps = current_nhaps;
            continue;
        }
        else if (previous_nhaps != current_nhaps) {
            HANDLE_ERROR(("line " + std::to_string(nloci_before_filter) +
                          " of " + filename + " has " + std::to_string(current_nhaps) +
                          " haps, but the previous line has " + std::to_string(previous_nhaps) + ".").c_str());
        }
        
        previous_nhaps = current_nhaps;
    }
    fin.clear(); // clear error flags
    fin.close();

    if(p.SKIP) { //prefilter all sites < MAF
        LOG(ARG_SKIP << " set. Removing all variants < " << p.MAF << ".");
        LOG("Removed "<<skiplist.size() <<" low frequency all variants.");

    }else{
        LOG(ARG_KEEP << " set. NOT removing variants < " << p.MAF << ".");
    }
    
    if(p.UNPHASED){
        if (current_nhaps % 2 != 0)
        {
            HANDLE_ERROR("Number of haplotypes must be even for unphased.");
        }
    }
    current_nhaps = (p.UNPHASED) ? current_nhaps/2 : current_nhaps;
    
    LOG("Loading " << current_nhaps << " haplotypes and " << nloci_before_filter - skiplist.size() << " loci...");
    hapData.initHapData(current_nhaps, nloci_before_filter-skiplist.size());
    hapData.skipQueue = queue<int>();

    //PHASE 2: Open input File 2nd time To Load into Data Structure
    fin.open(filename.c_str());
    if (fin.fail())
    {
        HANDLE_ERROR("Failed to open " + filename + " for reading.");
    }
    getline(fin, line);
    stringstream ss;
    char allele1;
    int locus_after_filter = 0;
    for (int locus = 0; locus < nloci_before_filter; locus++)
    {
        if(!skiplist.empty()){
            if(skiplist.front() == locus){
                skiplist.pop();
                hapData.skipQueue.push(locus);
                getline(fin, line);
                continue;
            }
        }

        ss.str(line);
        line.erase(std::remove_if(line.begin(), line.end(), [](char c) {
                return std::isspace(static_cast<unsigned char>(c));
            }), line.end());

        int cols = (p.UNPHASED) ? current_nhaps*2 : current_nhaps;

        if (p.UNPHASED){
            
            
            for (int hap = 0; hap < cols; hap++){ 
                if (hap % 2 == 0){
                    if(line[hap]=='1'  && line[hap+1]=='1'){
                        hapData.addAllele2(locus_after_filter, hap/2);
                    }else if ( (line[hap]=='1' && line[hap+1]=='0') || (line[hap]=='0' && line[hap+1]=='1') ){ // ==1
                        hapData.addAllele1(locus_after_filter, hap/2);
                    }
                }
            }
        }else{ // PHASED
            for (int hap = 0; hap < current_nhaps; hap++)
            {   
                ss >> allele1;
                if (allele1 != '0' && allele1 != '1')
                {
                    HANDLE_ERROR("Alleles must be coded 0/1 only.");
                }
                if(allele1=='1'){
                    hapData.addAllele1(locus_after_filter, hap);
                }
            }
        }

        locus_after_filter++;
        getline(fin, line);
    }
    fin.close();
}


//tested-> unphased - lowmem, phased - lowmem
void HapMap::readHapDataVCF(string filename, HapData& hapData)
{
    //RELATED FLAGS: ARG_SKIP, ARG_KEEP, ARG_UNPHASED, ARG_LOW_MEM, ARG_MAF, ARG_MISSING

    initParamsInHap(hapData);

    if(p.MISSING_ALLOWED){
        LOG("WARNING: Missing entries allowed: will impute missing entries if less than threshold.");
        readHapDataVCFMissing(filename, hapData);
        return;
    }
    
    igzstream fin;
    queue<int> skiplist;

    // PHASE 1: Counting so that inititalization is smooth
    LOG("Opening " << filename << " to count number of haplotypes and loci...");

    fin.open(filename.c_str());

    if (fin.fail())
    {
        HANDLE_ERROR("Failed to open " + filename + " for reading.");
    }

    int numMapCols = 9;
    string line;
    int nloci_before_filtering = 0;
    int prev_ngts = -1;
    int current_ngts = 0;

    int skipcount = 0;

    int num_meta_data_lines = 0;
    while (getline(fin, line))  //Counts number of haps (cols) and number of loci (rows)
    {
        if (line[0] == '#') {
            num_meta_data_lines++;
            continue;
        }
        nloci_before_filtering++;
        current_ngts = countFields(line) - numMapCols;

        /********/
        string junk;
        string chr;
        char allele1, allele2, separator;
        std::stringstream ss(line);
        int number_of_1s = 0;
        int number_of_2s = 0;

        for (int i = 0; i < numMapCols; i++) {
            ss >> junk;
            if(i==0){
                chr = junk;
            }
        }
        for (int field = 0; field < current_ngts; field++)
        {
            ss >> junk;
            allele1 = junk[0];
            separator = junk[1];
            allele2 = junk[2];

            if(!p.MISSING_ALLOWED){
                if ( (allele1 != '0' && allele1 != '1') || (allele2 != '0' && allele2 != '1') )
                {
                    HANDLE_ERROR("Alleles must be coded 0 or 1 only. Found alleles " << allele1 << separator << allele2);
                }
            }

            if(separator != '|' && !p.UNPHASED){
               HANDLE_ERROR("Unphased entries detected (| is used). Make sure you run with --unphased flag for correct results.");
            }

            if(p.UNPHASED){
                if (allele1 == '1' && allele2 == '1'){
                    number_of_2s++;
                }
                else if (allele1 == '1' || allele2 == '1'){
                    number_of_1s++;
                }
            }else{
                if(allele1 == '1'){
                    number_of_1s++;
                }
                if(allele2 == '1'){
                    number_of_1s++;
                }
            }
        }

        int nalleles_per_loc = current_ngts*2;
        bool skip_due_to_maf = shouldSkipLocus(number_of_1s, number_of_2s, nalleles_per_loc);
        
        // bool skipreason2 = false;
        // if(p.MULTI_CHR){ 
        //     if(chr_set.empty()){
        //         //cerr<<"WARNING: No chromosome list provided. Running analysis on all chromosomes.\n";
        //     }else{
        //         skipreason2 = (chr_set.find(chr) != chr_set.end());
        //     }
        // }

        if(!p.CALC_XP &&  !p.CALC_XPNSL){
            if ( skip_due_to_maf ) {
                skiplist.push(nloci_before_filtering-1);
                skipcount++;
            } 
        }
        
        /*********/
        if (prev_ngts < 0)
        {
            prev_ngts = current_ngts;
            continue;
        }
        else if (prev_ngts != current_ngts)
        {
            HANDLE_ERROR("line " + std::to_string(nloci_before_filtering) +
                          " of " + filename + " has " + std::to_string(current_ngts) +
                          " fields, but the previous line has " + std::to_string(prev_ngts) + ".");
        }
        prev_ngts = current_ngts;
    }

    fin.clear(); // clear error flags
    //fin.seekg(fileStart);
    fin.close();


    //PHASE 2: Load according to first pass information
    fin.open(filename.c_str());

    if (fin.fail())
    {
        HANDLE_ERROR("Failed to open " + filename + " for reading.");
    }

    if(p.SKIP && !p.CALC_XP &&  !p.CALC_XPNSL){
        LOG(ARG_SKIP << " set. Removing all variants < " << p.MAF << ".\n");
    }
   
    int nhaps = p.UNPHASED ? (current_ngts ) : (current_ngts ) * 2;
    LOG("Loading " << nhaps << " haplotypes and " << nloci_before_filtering-skipcount << " loci... skipped "<<skipcount << "loci \n");
    hapData.initHapData(nhaps, nloci_before_filtering-skipcount);

    string junk;
    char allele1, allele2, separator;
    bool skipLine = false; // to skip metadata lines

    hapData.skipQueue = skiplist; 
    int nloci_after_filtering = 0;

    for (int locus = 0; locus < nloci_before_filtering; locus++)
    {
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

        for (int field = 0; field <  current_ngts ; field++)
        {
            fin >> junk;
            allele1 = junk[0];
            separator = junk[1];

            allele2 = junk[2];

            if(p.UNPHASED){
                if (allele1 == '1' && allele2 == '1'){
                    hapData.addAllele2(nloci_after_filtering, field);
                }
                else if (allele1 == '1' || allele2 == '1'){
                    hapData.addAllele1(nloci_after_filtering, field);
                }
            }else{ // phased
                if(allele1 == '1'){
                    hapData.addAllele1(nloci_after_filtering, 2*field);
                }
                if(allele2 == '1'){
                    hapData.addAllele1(nloci_after_filtering, 2*field+1);
                }   
            }
        }
        nloci_after_filtering++;
    }
    

    if(p.SKIP){
        LOG("Removed " << skipcount << " low frequency variants from haplotype data.");
    }
    LOG("====");
    fin.close();

    // // // DEBUG
    // for(int locus_after_filter = 0; locus_after_filter < 1; locus_after_filter++){
    //     cout<<locus_after_filter<<"::: ";
    //     for(int i=0; i< hapEntries[locus_after_filter].positions.size(); i++){
    //         cout<<hapEntries[locus_after_filter].positions[i]<<" ";
    //     }
    //     cout<<endl;
        
    //     cout<<locus_after_filter<<":::*";
    //     for(int i=0; i< hapEntries[locus_after_filter].positions2.size(); i++){
    //         cout<<hapEntries[locus_after_filter].positions2[i]<<" ";
    //     }
    //     cout<<endl;
    // }
}




//tested-> unphased - lowmem, phased - lowmem
//RELATED FLAGS: ARG_SKIP, ARG_KEEP, ARG_UNPHASED, ARG_LOW_MEM, ARG_MAF, ARG_MISSING
void HapMap::readHapDataVCFMissing(string filename, HapData& hapData)
{
    initParamsInHap(hapData);
    float ALLOWED_MISSING_PERC = 1;
    
    igzstream fin;

    std::unique_ptr<std::vector<int> > vp1(new std::vector<int>());
    std::unique_ptr<std::vector<int> > vp2(new std::vector<int>());

    vector<int>& number_of_1s_per_loci = *vp1;
    vector<int>& number_of_2s_per_loci = *vp2;
    queue<int> skiplist;

    // PHASE 1: Counting so that inititalization is smooth
    LOG("Opening " << filename << " to count number of haplotypes and loci...");
    fin.open(filename.c_str());

    if (fin.fail())
    {
        HANDLE_ERROR("Failed to open " + filename + " for reading.");
    }

    int numMapCols = 9;
    string line;
    int nloci_before_filtering = 0;
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

        /********/
        string junk;
        char allele1, allele2, separator;
        std::stringstream ss(line);
        int number_of_1s = 0;
        int number_of_2s = 0;

        for (int i = 0; i < numMapCols; i++) {
            ss >> junk;
        }

        int missing_count = 0; // reset at every locus
        for (int field = 0; field < current_nhaps; field++)
        {
            ss >> junk;
            allele1 = junk[0];
            separator = junk[1];
            allele2 = junk[2];
            

            if(p.UNPHASED){
                char allele = '0';
                if (allele1 == '1' && allele2 == '1'){
                    number_of_2s++;
                }
                else if (allele1 == '1' || allele2 == '1'){
                    number_of_1s++;
                }
                if ( (allele1 != '0' && allele1 != '1') || (allele2 != '0' && allele2 != '1') )
                {
                    // if(allele1 == '.' || allele2 == '.'){
                    //     continue;
                    // }

                    LOG("WARNING: Alleles not coded 0/1, treating this as missing site.");


                    if(p.MISSING_MODE=="RANDOM"){
                        allele1 = randomZeroOrOne();
                        allele2 = randomZeroOrOne();
                    }else if(p.MISSING_MODE=="ONE_IMPUTE"){
                        allele1 = '1';
                        allele2 = '1';
                    }else if(p.MISSING_MODE=="ZERO_IMPUTE"){
                        allele1 = '0';
                        allele2 = '0';
                    }else if(p.MISSING_MODE=="NO_IMPUTE"){
                        missing_count+= 1;
                    }
                }
            }else{
                if( (allele1 != '0' && allele1 != '1') || (allele2 != '0' && allele2 != '1')){
                    if(p.MISSING_MODE=="RANDOM"){
                        allele1 = randomZeroOrOne();
                        allele2 = randomZeroOrOne();
                    }else if(p.MISSING_MODE=="ONE_IMPUTE"){
                        allele1 = '1';
                        allele2 = '1';
                    }else if(p.MISSING_MODE=="ZERO_IMPUTE"){
                        allele1 = '0';
                        allele2 = '0';
                        missing_count = 0;
                    }else if(p.MISSING_MODE=="NO_IMPUTE"){
                        missing_count+= 2;
                    }
                }
                if(allele1 == '1'){
                    number_of_1s++;
                }
                if(allele2 == '1'){
                    number_of_1s++;
                }
                
                if( (allele1 != '0' && allele1 != '1') || (allele2 != '0' && allele2 != '1')){
                    missing_count+= 2;
                }
            }
        }

        int derived_allele_count = (p.UNPHASED? (number_of_1s + number_of_2s*2) : number_of_1s);
      
        // --skip-low-freq filtering based on MAF
        //missing filter, being on the safe side
        if ( p.SKIP ){ // works for both phased and unphased
            int ancestral_allele_count = 0;
            ancestral_allele_count = current_nhaps*2 - derived_allele_count;
            bool derived_lower = (derived_allele_count+missing_count)*1.0/(current_nhaps*2) < p.MAF;
            bool ancestral_lower = (ancestral_allele_count+missing_count)*1.0/(current_nhaps*2) < p.MAF;
            bool missing_cutoff_crossed = 1.0*missing_count/(current_nhaps*2) > ALLOWED_MISSING_PERC;
            if (derived_lower || ancestral_lower || missing_cutoff_crossed) {
                //cout<<"Skipping cause Derived "<<derived_allele_count<<" Missing "<<missing_count<<" Total "<<current_nhaps*2<<endl;
                LOG("Skipping cause Derived " + std::to_string(derived_allele_count) + " Missing " + std::to_string(missing_count) + " Total " + std::to_string(current_nhaps*2));
                skiplist.push(nloci_before_filtering-1);
                skipcount++;
            }
        } 

        /*********/
        if (previous_nhaps < 0)
        {
            previous_nhaps = current_nhaps;
            continue;
        }
        else if (previous_nhaps != current_nhaps)
        {
            HANDLE_ERROR("line " + std::to_string(nloci_before_filtering) +
                          " of " + filename + " has " + std::to_string(current_nhaps) +
                          " fields, but the previous line has " + std::to_string(previous_nhaps) + ".");
        }
        previous_nhaps = current_nhaps;
    }

    fin.clear(); // clear error flags
    //fin.seekg(fileStart);
    fin.close();


    //PHASE 2: Load according to first pass information
    fin.open(filename.c_str());

    if (fin.fail())
    {
        HANDLE_ERROR("Failed to open " + filename + " for reading.");
    }

    if(p.SKIP){
        LOG(ARG_SKIP << " set. Removing all variants < " << p.MAF << ".");
    }
   
    int nhaps = p.UNPHASED ? (current_nhaps ) : (current_nhaps ) * 2;
    LOG("Loading " << nhaps << " haplotypes and " << nloci_before_filtering-skipcount << " loci... skipped "<<skipcount << "loci \n");
    hapData.initHapData(nhaps, nloci_before_filtering-skipcount);

    string junk;
    char allele1, allele2, separator;
    bool skipLine = false; // to skip metadata lines

    hapData.skipQueue = skiplist; 
    int nloci_after_filtering = 0;

    for (int locus = 0; locus < nloci_before_filtering; locus++)
    {
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
        

        for (int field = 0; field <  current_nhaps ; field++)
        {
            fin >> junk;
            allele1 = junk[0];
            separator = junk[1];

            allele2 = junk[2];
            if ( (allele1 != '0' && allele1 != '1') || (allele2 != '0' && allele2 != '1') )
            {

                LOG("WARNING: Alleles must be coded 0/1 only. Found "<< allele1 << " " << allele2);
                allele1 = '?';
                allele2 = '?';
            }

            if(p.UNPHASED){
                char allele = '0';
                if (allele1 == '1' && allele2 == '1'){
                    hapData.addAllele2(nloci_after_filtering, field);
                }
                else if ((allele1 == '1' && allele2 == '0') || (allele2 == '1' && allele1 == '0')){
                    hapData.addAllele1(nloci_after_filtering, field);
                }else if(allele1 == '?' || allele2 == '?'){
                    if( (allele1 != '0' && allele1 != '1') || (allele2 != '0' && allele2 != '1')){
                        if(p.MISSING_MODE=="RANDOM"){
                            allele1 = randomZeroOrOne();
                            allele2 = randomZeroOrOne();
                        }else if(p.MISSING_MODE=="ONE_IMPUTE"){
                            allele1 = '1';
                            allele2 = '1';
                        }else if(p.MISSING_MODE=="ZERO_IMPUTE"){
                            allele1 = '0';
                            allele2 = '0';
                        }else if(p.MISSING_MODE=="NO_IMPUTE"){
                            hapData.addAlleleMissing(nloci_after_filtering, field);
                        }
                    }
                }
            }else{ // phased
                if(allele1 == '1'){
                    hapData.addAllele1(nloci_after_filtering, 2*field);
                }
                if(allele2 == '1'){
                    hapData.addAllele1(nloci_after_filtering, 2*field + 1);
                }   
                if(allele1 == '?' || allele2 == '?'){
                    if( (allele1 != '0' && allele1 != '1') || (allele2 != '0' && allele2 != '1')){
                        if(p.MISSING_MODE=="RANDOM"){
                            allele1 = randomZeroOrOne();
                            allele2 = randomZeroOrOne();
                        }else if(p.MISSING_MODE=="ONE_IMPUTE"){
                            allele1 = '1';
                            allele2 = '1';
                        }else if(p.MISSING_MODE=="ZERO_IMPUTE"){
                            allele1 = '0';
                            allele2 = '0';
                        }else if(p.MISSING_MODE=="NO_IMPUTE"){
                            hapData.addAlleleMissing(nloci_after_filtering, 2*field);
                            hapData.addAlleleMissing(nloci_after_filtering, 2*field+1);
                            // hapData.hapEntries[nloci_after_filtering].xorbitset->set_bit(2 * field);
                            // hapData.hapEntries[nloci_after_filtering].xorbitset->set_bit(2 * field + 1);
                            // hapData.hapEntries[nloci_after_filtering].xorbitset->num_1s+= 2;
                        }
                    }
                }
            }
        }
        nloci_after_filtering++;
    }

    if(p.SKIP){
        LOG("Removed " << skipcount << " low frequency variants from haplotype data.");
    }

    fin.close();

    // // // DEBUG
    // for(int locus_after_filter = 0; locus_after_filter < 1; locus_after_filter++){
    //     cout<<locus_after_filter<<"::: ";
    //     for(int i=0; i< hapEntries[locus_after_filter].positions.size(); i++){
    //         cout<<hapEntries[locus_after_filter].positions[i]<<" ";
    //     }
    //     cout<<endl;
        
    //     cout<<locus_after_filter<<":::*";
    //     for(int i=0; i< hapEntries[locus_after_filter].positions2.size(); i++){
    //         cout<<hapEntries[locus_after_filter].positions2[i]<<" ";
    //     }
    //     cout<<endl;
    // }
}

// Decide whether to skip a locus based on MAF and missing data
bool HapMap::shouldSkipLocus(int number_of_1s, int number_of_2s, int nalleles_per_loc, int missing_count)
{
    // Compute the count of derived alleles depending on phased/unphased data
    int derived_allele_count = (p.UNPHASED
        ? (number_of_1s + number_of_2s * 2)  // unphased: 1s = het, 2s = homozygous derived
        : number_of_1s);                    // phased: count of 1s only

    // Flags for filtering
    bool skip_due_to_maf = false;     // Skip based on minor allele frequency threshold
    
    // Check for low-frequency alleles (MAF filtering)
    if (p.SKIP) {
        double derived_freq = derived_allele_count * 1.0 / nalleles_per_loc;
        double ancestral_freq = 1.0 - derived_freq;

        skip_due_to_maf = (derived_freq < MIN_MAF_CUTOFF || ancestral_freq < MIN_MAF_CUTOFF);
    }
    
    return skip_due_to_maf;

    /* UNCOMMENT if missing data is allowed
    bool skip_due_to_missing = false; // Skip based on missingness threshold
    if (p.MISSING_ALLOWED) {
        const double ALLOWED_MISSING_PERC = 1.0; // 100% allowed unless explicitly used
        int ancestral_allele_count = nalleles_per_loc - derived_allele_count;

        // Check if MAF threshold would be violated when accounting for missing data
        bool derived_too_rare = (derived_allele_count + missing_count) * 1.0 / nalleles_per_loc < MIN_MAF_CUTOFF;
        bool ancestral_too_rare = (ancestral_allele_count + missing_count) * 1.0 / nalleles_per_loc < MIN_MAF_CUTOFF;

        if (p.SKIP) {
            skip_due_to_missing = derived_too_rare || ancestral_too_rare;
        }

        // Check if missingness exceeds allowed threshold
        bool too_much_missing = (1.0 * missing_count / nalleles_per_loc) > ALLOWED_MISSING_PERC;
        skip_due_to_missing = skip_due_to_missing || too_much_missing;
    }
    return skip_due_to_maf || skip_due_to_missing;
    */


}


// bool HapMap::skipThisLocus(int number_of_1s, int number_of_2s, int nalleles_per_loc, int missing_count){ 
//     //nalleles same for phased and unphased, tot n of columns

//     int derived_allele_count = (p.UNPHASED? (number_of_1s + number_of_2s*2) : number_of_1s);
    
//     // if ((derived_allele_count)*1.0/(nalleles_per_loc) == 0){
//     //     return true;
//     // }

//     // if( 1-(derived_allele_count*1.0/(nalleles_per_loc)) == 0 ){
//     //     return true;
//     // }

//     bool skipreason1 = false; // --skip-low-freq filtering based on MAF
//     bool skipreason2 = false; // missing

//     double ALLOWED_MISSING_PERC = 1;
//     if(p.MISSING_ALLOWED){
//         int ancestral_allele_count = nalleles_per_loc - derived_allele_count;
//         bool derived_lower = (derived_allele_count+missing_count)*1.0/(nalleles_per_loc) < this->MIN_MAF_CUTOFF ;
//         bool ancestral_lower = (ancestral_allele_count+missing_count)*1.0/(nalleles_per_loc) < this->MIN_MAF_CUTOFF ;
//         if(p.SKIP){
//                 skipreason2 = (derived_lower || ancestral_lower);
//         }

//         bool missing_cutoff_crossed = 1.0*missing_count/nalleles_per_loc > ALLOWED_MISSING_PERC;
//         skipreason2 = (skipreason2 || missing_cutoff_crossed);
//         // if (missing_cutoff_crossed) {
//         //     //cout<<"Skipping cause Derived "<<derived_allele_count<<" Missing "<<missing_count<<" Total "<<current_nhaps*2<<endl;
//         //     skipreason2 = true;
//         // }
//     }
    
//     //if(!p.MISSING_ALLOWED){
//         skipreason1 = (p.SKIP && (derived_allele_count*1.0/(nalleles_per_loc) < this->MIN_MAF_CUTOFF || 1-(derived_allele_count*1.0/(nalleles_per_loc)) < this->MIN_MAF_CUTOFF ));    
//     //}
//     return (skipreason1 ||  skipreason2 );
// }





// TPED format description
// 1. Chromosome (string)
// 2. SNP ID (rs number: string)
// 3. Genetic distance (cM: float)
// 4. Physical position (int)
// 5. Allele 1: Sample 1 (0 or 1: char)
// 6. Allele 2: Sample 1 (0 or 1: char)
// 7. Allele 1: Sample 2 (0 or 1: char)
// 8. Allele 2: Sample 2 (0 or 1: char)
// ...
void HapMap::readHapDataTPED(string filename, HapData &hapData)
{
    initParamsInHap(hapData);
    igzstream fin;
    
    // PHASE1: opening file first time and counting number of loci and haps, and determine which loci to skip
    LOG("Opening " << filename << " to count number of haplotypes and loci...");

    fin.open(filename.c_str());
    if (fin.fail())
    {
        HANDLE_ERROR("Failed to open " + filename + " for reading.");
    }

    int numMapCols = 4;
    //int fileStart = fin.tellg();
    string line;
    int nloci_before_filter = 0; // number of loci before filtering (filtering happens when SKIP is enabled)
    int previous_nhaps = -1;
    int current_nhaps = 0;
    
    queue<int> skiplist;

    //Counts number of haps (cols) and number of loci (rows)
    //if any lines differ, send an error message and throw an exception
    while (getline(fin, line))
    {
        //getline(fin,line);
        //if(fin.eof()) break;
        nloci_before_filter++; // each line contributes to a locus
        current_nhaps = countFields(line) - numMapCols; // here, current_nhaps (count of fields) is the number of haps 
        
        stringstream ss;
        ss.str(line);

        char allele1, allele2;
        string junk;
        
        for(int i = 0; i < numMapCols; i++){
            ss >> junk;
        }

        std::ostringstream rest;
        while (ss >> junk) {
            rest << junk ;
        }

        line.clear();
        string line_haps = rest.str();
        line_haps.erase(std::remove_if(line_haps.begin(), line_haps.end(), [](char c) {
            return std::isspace(static_cast<unsigned char>(c));
        }), line_haps.end());  // erase spaces inside and outside


        int number_of_2s = 0; // for unphased
        int number_of_1s = 0;
        for (int hap = 0; hap < current_nhaps; hap++){ 
            if(p.UNPHASED){
                if (hap % 2 == 0){
                    if(line_haps[hap]=='1'  && line_haps[hap+1]=='1'){
                        number_of_2s++;
                    }else if ( (line_haps[hap]=='1' && line_haps[hap+1]=='0') || (line_haps[hap]=='0' && line_haps[hap+1]=='1') ){ // ==1
                        number_of_1s++;
                    }
                }
            }else{
                if(line_haps[hap]=='1'){
                    number_of_1s++;
                }
            }
        }
        if(shouldSkipLocus(number_of_1s, number_of_2s, current_nhaps)){
            skiplist.push(nloci_before_filter-1);
        }

        if (previous_nhaps < 0)
        {
            previous_nhaps = current_nhaps;
            continue;
        }
        else if (previous_nhaps != current_nhaps) {
            HANDLE_ERROR(("line " + std::to_string(nloci_before_filter) +
                          " of " + filename + " has " + std::to_string(current_nhaps) +
                          " haps, but the previous line has " + std::to_string(previous_nhaps) + ".").c_str());
        }
        
        previous_nhaps = current_nhaps;
    }

    fin.clear(); // clear error flags
    //fin.seekg(fileStart);
    fin.close();


    if(p.SKIP) { //prefilter all sites < MAF
        LOG(ARG_SKIP << " set. Removing all variants < " << p.MAF << ".");
        LOG("Removed "<<skiplist.size() <<" low frequency all variants.");

    }else{
        LOG(ARG_KEEP << " set. NOT removing variants < " << p.MAF << ".");
    }
    
    if(p.UNPHASED){
        if (current_nhaps % 2 != 0)
        {
            HANDLE_ERROR("Number of haplotypes must be even for unphased.");
        }
    }
    current_nhaps = (p.UNPHASED) ? current_nhaps/2 : current_nhaps;
    
    LOG("Loading " << current_nhaps << " haplotypes and " << nloci_before_filter - skiplist.size() << " loci...");
    hapData.initHapData(current_nhaps, nloci_before_filter-skiplist.size());
    hapData.skipQueue = queue<int>();

    
    // PHASE 2: Open again and load the data to HapData: now that we gathered some idea about the format of the file
    fin.open(filename.c_str());

    if (fin.fail())
    {
        HANDLE_ERROR("Failed to open " + filename + " for reading.");
    }

    
    string junk;
    int locus_after_filter = 0;
    for (int locus = 0; locus < nloci_before_filter; locus++)
    {
        getline(fin, line);
        if(!skiplist.empty()){
            if(skiplist.front() == locus){
                skiplist.pop();    
                hapData.skipQueue.push(locus); // make a copy to skipQueue
                continue;
            }
        }

        istringstream ss(line);
        for (int i = 0; i < numMapCols; i++)
        {
            ss >> junk; // skip first 4 columns
        }

        std::ostringstream rest;

        int cols = current_nhaps; 
        
        if(p.UNPHASED) {
            cols *= 2;
        }
        for (int i = 0; i < cols; i++)
        {
            ss >> junk; 
            rest << junk << "";
        }

        string line_haps = rest.str();
        
        line_haps.erase(std::remove_if(line_haps.begin(), line_haps.end(), [](char c) {
            return std::isspace(static_cast<unsigned char>(c));
        }), line_haps.end());  // erase spaces inside and outside

        if (p.UNPHASED){
            for (int hap = 0; hap < cols; hap++){ 
                if (hap % 2 == 0){
                    if(line_haps[hap]=='1'  && line_haps[hap+1]=='1'){
                        hapData.addAllele2(locus_after_filter, hap/2);
                    }else if ( (line_haps[hap]=='1' && line_haps[hap+1]=='0') || (line_haps[hap]=='0' && line_haps[hap+1]=='1') ){ // ==1
                        hapData.addAllele1(locus_after_filter, hap/2);
                    }
                }
            }
        }else{ // PHASED
            for (int hap = 0; hap < current_nhaps; hap++)
            {   
                char allele1 = line_haps[hap];
                if (allele1 != '0' && allele1 != '1')
                {
                    HANDLE_ERROR("Alleles must be coded 0/1 only.");
                }
                if(allele1=='1'){
                    hapData.addAllele1(locus_after_filter, hap);
                }
            }
        }
        locus_after_filter++;
        line.clear();
        ss.clear();
        rest.clear();
    }

    fin.close();
}

void HapMap::loadHapMapData(){
        
    mapData = std::make_unique<MapData>(); 
    hapData = std::make_unique<HapData>();

    // PHASE 2:  Load data to HapData and MapData
    if (p.CALC_XP || p.CALC_XPNSL){
        hapData2 = std::make_unique<HapData>();
    }   

    auto start_reading = std::chrono::high_resolution_clock::now();

    if (p.TPED)
    {
        readHapDataTPED(p.tpedFilename, *hapData);
        if (p.CALC_XP || p.CALC_XPNSL)
        {
            readHapDataTPED(p.tpedFilename, *hapData2);
            if (hapData->nloci != hapData2->nloci)
            {
                HANDLE_ERROR("Haplotypes from " + p.tpedFilename + " and " + p.tpedFilename2 + " do not have the same number of loci.");
            }
        }
        int exp_haps = (p.UNPHASED)? hapData->nhaps*2 : hapData->nhaps;
        mapData->readMapDataTPED(p.tpedFilename, hapData->nloci, exp_haps, p.USE_PMAP, hapData->skipQueue);
    }
    else if (p.VCF) {
        if (p.CALC_XP || p.CALC_XPNSL)
        {
            readHapDataVCF(p.vcfFilename, *hapData);
            readHapDataVCF(p.vcfFilename2, *hapData2);
            if (hapData->nloci != hapData2->nloci)
            {
                HANDLE_ERROR("Haplotypes from " + p.vcfFilename + " and " + p.vcfFilename2 + " do not have the same number of loci.");
            }
        }else{
            readHapDataVCF(p.vcfFilename, *hapData);
        }
        if(!p.CALC_NSL && !p.CALC_XPNSL && !p.USE_PMAP) {
            //int exp_haps = (p.UNPHASED)? hapData->nhaps*2 : hapData->nhaps;
            mapData->readMapData(p.mapFilename, hapData->nloci, p.USE_PMAP, hapData->skipQueue);
        }
        else{//Load physical positions
            mapData->readMapDataVCF(p.vcfFilename, hapData->nloci, hapData->skipQueue);
        }
    }else if(p.THAP){
        readHapDataTHAP(p.thapFilename, *hapData);
        if (p.CALC_XP || p.CALC_XPNSL)
        {
            readHapDataTHAP(p.thapFilename2, *hapData2);
            if (hapData->nloci != hapData2->nloci)
            {
                HANDLE_ERROR("Haplotypes from " + p.thapFilename + " and " + p.thapFilename2 + " do not have the same number of loci.");
            }
        }
        mapData->readMapData(p.mapFilename, hapData->nloci, p.USE_PMAP, hapData->skipQueue);
    }else
    {
        readHapDataMSHap(p.hapFilename, *hapData);
        if (p.CALC_XP || p.CALC_XPNSL)
        {
            readHapDataMSHap(p.hapFilename2, *hapData2);
            if (hapData->nloci != hapData2->nloci)
            {
                HANDLE_ERROR("Haplotypes from " + p.hapFilename + " and " + p.hapFilename2 + " do not have the same number of loci.");
            }
        }
        mapData->readMapData(p.mapFilename, hapData->nloci, p.USE_PMAP, hapData->skipQueue);
    }

    // PHASE 2:  XOR
    // if(p.benchmark_flag=="XOR")
    //     hapData->xor_for_phased_and_unphased();

    auto end_reading = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> read_duration =  end_reading - start_reading;
    LOG("Input file loaded in "+to_string(read_duration.count())+" s.");
    LOG("=====");

    // DEBUG::: mapData.print();

    // if(mapData->chr_list.size() > 1 && p.MULTI_CHR == false){
    //     cerr<<"ERROR: Input file contains multiple chromosomes, but no chromosome list is provided. Use --multi-chr flag to specify chromosomes.\n";
    // }

    // Check if map is in order
    for (int i = 1; i < mapData->nloci; i++) {
        if ( mapData->mapEntries[i].physicalPos <  mapData->mapEntries[i-1].physicalPos ) {
            string ERROR_MSG = "Variant physical position must be monotonically increasing.\n";
            ERROR_MSG += "\t" + mapData->mapEntries[i].locusName + " " + std::to_string(mapData->mapEntries[i].physicalPos) + " appears after";
            ERROR_MSG += "\t" + mapData->mapEntries[i-1].locusName + " " + std::to_string(mapData->mapEntries[i-1].physicalPos) + ".";
            HANDLE_ERROR(ERROR_MSG);
        }
        if ( !p.CALC_NSL && mapData->mapEntries[i].geneticPos  < mapData->mapEntries[i-1].geneticPos  ) {
            string ERROR_MSG = "Variant genetic position must be monotonically increasing.\n";
            ERROR_MSG += "\t" + mapData->mapEntries[i].locusName + " " + std::to_string(mapData->mapEntries[i].geneticPos) + " appears after";
            ERROR_MSG += "\t" + mapData->mapEntries[i-1].locusName + " " + std::to_string(mapData->mapEntries[i-1].geneticPos) + ".";
            HANDLE_ERROR(ERROR_MSG);
        }
    }
    
}
