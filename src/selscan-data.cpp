
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

    if(p.CALC_IHS || p.CALC_NSL){
        if(p.MULTI_MAF){
            for (const auto& param : ps) {
                if(param.MAF < MIN_MAF_CUTOFF){
                    this->MIN_MAF_CUTOFF = param.MAF;
                }
            }
        }
    }else{
        // if we are not calculating IHS or NSL, we can set the MAF cutoff to 0
        this->MIN_MAF_CUTOFF = 0;
    }
    // if IHS OR NSL othwerwise MINMAF==0
}




// ms hap format : number of rows = number of haplotypes, number of columns = number of loci
void HapMap::readHapDataMSHap(string filename, HapData &hapData)
{
    //HANDLE_ERROR("Bug found in ms-style hap format. This function is undergoing development.");
    initParamsInHap(hapData);

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
            
            if (skip_idx < static_cast<int>(skiplist.size())) {
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
    hapData.skipQueue = queue<size_t>();

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



VCFPass1Result HapMap::readHapDataVCF_pass1(string filename)
{
    VCFPass1Result result;

    igzstream fin;
    fin.open(filename.c_str());
    if (fin.fail()) {
        HANDLE_ERROR("Failed to open " + filename + " for reading.");
    }

    constexpr int numMapCols = 9;
    string line;

    int physpos = -1;
    int physpos_first_duplicated_id = -1;
    int skip_due_to_duplicate_pos = 0;

    auto find_field_end = [](const char* p, const char* end) -> const char* {
        while (p < end && *p != '\t') ++p;
        return p;
    };

    auto parse_int_fast = [](const char* s, const char* e) -> int {
        int x = 0;
        for (const char* p = s; p < e; ++p) {
            x = x * 10 + (*p - '0');
        }
        return x;
    };

    auto count_samples_from_header = [&](const string& header_line) -> int {
        int tabs = 0;
        for (char c : header_line)
            if (c == '\t') ++tabs;
        return (tabs + 1) - numMapCols;
    };

    auto is_snp_base = [](const char* s, const char* e) -> bool {
        if (e - s != 1) return false;
        char c = *s;
        return c=='A'||c=='C'||c=='G'||c=='T';
    };

    while (getline(fin, line))
    {
        if (line.empty()) continue;

        if (line[0] == '#') {
            if (line.rfind("#CHROM",0) == 0) {
                result.current_ngts = count_samples_from_header(line);
            }
            continue;
        }

        result.nloci_before_filtering++;

        // bool ALLOW_XP_LOCI_MISMATCH = false; // for now, assume all loci mismatch in XP mode, and rely on physpos-based skipping to identify shared loci. This is because in XP mode we want to keep all loci that are present in either VCF, even if they are not shared. So we will disable the NUM_LOCI_MISMATCH check and just rely on physpos-based skipping to identify shared loci.
        // if((p.CALC_XP || p.CALC_XPNSL) && ALLOW_XP_LOCI_MISMATCH ){
        //     result.physpos.push_back(physpos);
        // }

        const char* pcur = line.c_str();
        const char* end  = pcur + line.size();

        const char* field_starts[numMapCols];
        const char* field_ends[numMapCols];

        for (int col = 0; col < numMapCols; ++col) {
            field_starts[col] = pcur;
            field_ends[col]   = find_field_end(pcur,end);
            pcur = field_ends[col];
            if (pcur < end) ++pcur;
        }

        int new_physpos = parse_int_fast(field_starts[1], field_ends[1]);

        if (physpos == new_physpos)
            skip_due_to_duplicate_pos++;
        else {
            skip_due_to_duplicate_pos = 0;
            physpos_first_duplicated_id = result.nloci_before_filtering - 1;
        }

        physpos = new_physpos;

        bool skip_due_to_missing = false;
        bool skip_due_to_multiallelic = false;

        int number_of_1s = 0;
        int number_of_2s = 0;

        if (p.MULTI_ALLELIC) {
            if (!is_snp_base(field_starts[3],field_ends[3]) ||
                !is_snp_base(field_starts[4],field_ends[4])) {
                skip_due_to_multiallelic = true;
            }
        }

        int seen_ngts = 0;

        if (p.UNPHASED)
        {
            while (pcur < end && seen_ngts < result.current_ngts)
            {
                const char* g = pcur;

                if (end-g >= 3) {
                    char a1 = g[0];
                    char a2 = g[2];

                    int v1 = (a1=='1');
                    int v2 = (a2=='1');

                    int sum = v1+v2;

                    if (sum==2) number_of_2s++;
                    else if (sum==1) number_of_1s++;
                }

                while (pcur < end && *pcur!='\t') ++pcur;
                if (pcur < end) ++pcur;

                seen_ngts++;
            }
        }
        else
        {
            while (pcur < end && seen_ngts < result.current_ngts)
            {
                const char* g = pcur;

                if (end-g >= 3) {
                    char a1 = g[0];
                    char a2 = g[2];

                    number_of_1s += (a1=='1');
                    number_of_1s += (a2=='1');
                }

                while (pcur < end && *pcur!='\t') ++pcur;
                if (pcur < end) ++pcur;

                seen_ngts++;
            }
        }

        // -------- store counts per locus --------
        result.num_1s.push_back(number_of_1s);
        result.num_2s.push_back(number_of_2s);

        int nalleles_per_loc = result.current_ngts * 2;
        bool skip_due_to_maf = shouldSkipLocus(number_of_1s, number_of_2s, nalleles_per_loc);

        if (p.CALC_XP || p.CALC_XPNSL)
            skip_due_to_maf = false;

        if (skip_due_to_duplicate_pos == 1) {
            if (result.skiplist.empty() ||
                result.skiplist.back() != (size_t)physpos_first_duplicated_id)
            {
                result.skiplist.push(result.nloci_before_filtering-2);
                result.skipcount++;
            }
        }

        if (skip_due_to_maf || skip_due_to_multiallelic ||
            skip_due_to_missing || skip_due_to_duplicate_pos)
        {
            result.skiplist.push(result.nloci_before_filtering-1);
            result.skipcount++;
        }
    }

    fin.close();
    return result;
}


void HapMap::readHapDataVCF_pass2(string filename,  HapData& hapData, const VCFPass1Result& pass1)
{
    igzstream fin;
    fin.open(filename.c_str());
    if (fin.fail()) {
        HANDLE_ERROR("Failed to open " + filename + " for reading.");
    }

    constexpr int numMapCols = 9;   // fixed VCF columns before genotype fields
    string line;

    // Skip current field and land on first char of next field
    auto skip_to_next_field = [](const char* p, const char* end) -> const char* {
        while (p < end && *p != '\t') ++p;
        if (p < end) ++p;
        return p;
    };

    if (p.SKIP && !p.CALC_XP && !p.CALC_XPNSL) {
        LOG(ARG_SKIP << " set. Removing all variants < " << p.MAF << ".");
    }


    const int nhaps = p.UNPHASED ? pass1.current_ngts : pass1.current_ngts * 2;

    LOG("Loading " << nhaps
                   << " haplotypes and " << (pass1.nloci_before_filtering - pass1.skipcount)
                   << " loci."); 

    if(p.CALC_XP || p.CALC_XPNSL){
        LOG("XP mode: MAF-based filtering is disabled.");
        if(pass1.skipcount > 0)
            LOG(pass1.skipcount << " loci will be skipped.");

    }else{
        if (p.SKIP) {
            LOG(ARG_SKIP << " set. Removing all variants < " << p.MAF << ".");
            LOG("Removed " << pass1.skipcount << " low frequency variants from haplotype data.");
        } else {
            LOG(ARG_KEEP << " set. NOT removing variants < " << p.MAF << ".");
        }
    }

    // allocate exact size after filtering
    hapData.initHapData(nhaps, pass1.nloci_before_filtering - pass1.skipcount);

    // preserve your old behavior
    hapData.skipQueue = pass1.skiplist;

    // local copy because we consume it during pass 2
    queue<size_t> skiplist = pass1.skiplist;

    size_t locus = 0;                  // locus index in original file order
    size_t nloci_after_filtering = 0;  // locus index in filtered HapData

    while (getline(fin, line))
    {
        if (line.empty()) continue;
        if (line[0] == '#') continue;

        // skip loci identified in pass 1
        if (!skiplist.empty() && skiplist.front() == locus) {
            skiplist.pop();
            locus++;
            continue;
        }

        const char* pcur = line.c_str();
        const char* end  = pcur + line.size();

        // skip first 9 fixed VCF columns and land on first genotype field
        for (int col = 0; col < numMapCols; ++col) {
            pcur = skip_to_next_field(pcur, end);
        }

        if (p.UNPHASED)
        {
            // in unphased mode, each sample contributes one diploid state
            for (int field = 0; field < pass1.current_ngts; ++field)
            {
                // assume GT starts at first 3 chars: 0|1, 1/0, 0|1:...
                const int a1 = (pcur[0] == '1');
                const int a2 = (pcur[2] == '1');
                const int sum = a1 + a2;

                if (sum == 2) {
                    hapData.addAllele2(nloci_after_filtering, field);
                } else if (sum == 1) {
                    hapData.addAllele1(nloci_after_filtering, field);
                }

                while (pcur < end && *pcur != '\t') ++pcur;
                if (pcur < end) ++pcur;
            }
        }
        else
        {
            // in phased mode, each sample contributes two haplotypes
            for (int field = 0; field < pass1.current_ngts; ++field)
            {
                const int a1 = (pcur[0] == '1');
                const int a2 = (pcur[2] == '1');

                if (a1) hapData.addAllele1(nloci_after_filtering, 2 * field);
                if (a2) hapData.addAllele1(nloci_after_filtering, 2 * field + 1);

                while (pcur < end && *pcur != '\t') ++pcur;
                if (pcur < end) ++pcur;
            }
        }

        nloci_after_filtering++;
        locus++;
    }

    if(p.CALC_XP || p.CALC_XPNSL){
        LOG("Done removing " << pass1.skipcount << " loci due to duplicate positions.");
    }else{
        if (p.SKIP) {
            LOG("Done removing " << pass1.skipcount << " low frequency variants from haplotype data.");
        }
    }
    

    LOG("====");
    fin.close();


}


void HapMap::readHapDataVCFXP(string filename, string filename2, HapData& hapData, HapData& hapData2)
{
    // RELATED FLAGS:
    // ARG_SKIP, ARG_KEEP, ARG_UNPHASED, ARG_LOW_MEM, ARG_MAF, ARG_MISSING

    initParamsInHap(hapData);
    initParamsInHap(hapData2);


    LOG("Opening " << filename << " to count number of haplotypes and loci...");

    VCFPass1Result pass1_h1 = readHapDataVCF_pass1(filename);
    VCFPass1Result pass1_h2 = readHapDataVCF_pass1(filename2);
    // if(!pass1_h1.skiplist.empty() || !pass1_h2.skiplist.empty()){
    //     HANDLE_ERROR("XP should have no skipped loci due to MAF filtering.");
    // }
    //but skip maybe due to duplicate positions

    bool NUM_LOCI_MISMATCH = (pass1_h1.nloci_before_filtering != pass1_h2.nloci_before_filtering);
    // for (int pid = 0; pid < pass1_h1.nloci_before_filtering ; pid++) {
    //     if (pass1_h1.physpos[pid] != pass1_h2.physpos[pid]) {
    //         NUM_LOCI_MISMATCH = true;
    //         break;
    //     }
    // }
    if(NUM_LOCI_MISMATCH) {
        LOG("ERROR: The two VCF files have different sets of loci.");
        cout<<"File 1: " << filename << " has " << pass1_h1.nloci_before_filtering << " loci, with " << pass1_h1.skipcount << " skipped due to MAF or duplicates.\n";
        cout<<"File 2: " << filename2 << " has " << pass1_h2.nloci_before_filtering << " loci, with " << pass1_h2.skipcount << " skipped due to MAF or duplicates.\n";
       exit(EXIT_FAILURE);

//        LOG("WARNING: The two VCF files have different sets of loci. Will identify shared loci and skip non-shared ones.");
    } else {
        LOG("The two VCF files contain the same number of loci. Calculations are performed assuming they represent the same set of loci.");
    }
    
    // if (NUM_LOCI_MISMATCH) {            
    //     size_t i = 0;
    //     size_t j = 0;

    //     while (i < pass1_h1.physpos.size() && j < pass1_h2.physpos.size()) {

    //         if (pass1_h1.physpos[i] == pass1_h2.physpos[j]) {
    //             // value appears in both
    //             i++;
    //             j++;
    //         }
    //         else if (pass1_h1.physpos[i] < pass1_h2.physpos[j]) {
    //             // h1 value not present in h2
    //             pass1_h1.skiplist.push(i);
    //             i++;
    //         }
    //         else {
    //             // h2 value not present in h1
    //             pass1_h2.skiplist.push(j);
    //             j++;
    //         }
    //     }

    //     // remaining elements in h1
    //     while (i < pass1_h1.physpos.size()) {
    //         pass1_h1.skiplist.push(i);
    //         i++;
    //     }

    //     // remaining elements in h2
    //     while (j < pass1_h2.physpos.size()) {
    //         pass1_h2.skiplist.push(j);
    //         j++;
    //     }    
    // }

    //if the two vcfs differ in number of loci, we should handle that case. For now we assume they are the same.

    if (NUM_LOCI_MISMATCH)
    {
        vector<size_t> new_skip_h1;
        vector<size_t> new_skip_h2;


        const int total_alleles_h1 = 2 * pass1_h1.current_ngts; // a genotype is like 0|1, so 2 alleles per genotype
        const int total_alleles_h2 = 2 * pass1_h2.current_ngts;
        const int total_alleles_joint = total_alleles_h1 + total_alleles_h2;
        size_t i = 0;
        size_t j = 0;
       
        while (i < pass1_h1.physpos.size() || j < pass1_h2.physpos.size())
        {

            // Decide next genomic position to process
            int pos;
            if (i < pass1_h1.physpos.size() && j < pass1_h2.physpos.size())
            {
                pos = std::min(pass1_h1.physpos[i], pass1_h2.physpos[j]);
            }
            else if (i < pass1_h1.physpos.size())
            {
                pos = pass1_h1.physpos[i];
            }
            else
            {
                pos = pass1_h2.physpos[j];
            }

            // Find full run of this position in h1
            size_t i_start = i;
            while (i < pass1_h1.physpos.size() && pass1_h1.physpos[i] == pos)
            {
                i++;
            }
            size_t i_end = i;
            size_t count_h1 = i_end - i_start;

            // Find full run of this position in h2
            size_t j_start = j;
            while (j < pass1_h2.physpos.size() && pass1_h2.physpos[j] == pos)
            {
                j++;
            }
            size_t j_end = j;
            size_t count_h2 = j_end - j_start;

            // ----------------------------------------------------
            // Case 1: duplicated anywhere -> skip from everywhere
            // ----------------------------------------------------
            if (count_h1 > 1 || count_h2 > 1)
            {
                for (size_t k = i_start; k < i_end; k++)
                {
                    new_skip_h1.push_back(k);
                }
                for (size_t k = j_start; k < j_end; k++)
                {
                    new_skip_h2.push_back(k);
                }
                continue;
            }

            // ----------------------------------------------------
            // Case 2: exists only in h1
            // ----------------------------------------------------
            if (count_h1 == 1 && count_h2 == 0)
            {
                new_skip_h1.push_back(i_start);
                continue;
            }

            // ----------------------------------------------------
            // Case 3: exists only in h2
            // ----------------------------------------------------
            if (count_h1 == 0 && count_h2 == 1)
            {
                new_skip_h2.push_back(j_start);
                continue;
            }

            // ----------------------------------------------------
            // Case 4: shared exactly once in both files
            // ----------------------------------------------------
            if (count_h1 == 1 && count_h2 == 1)
            {

                // int derived_h1, derived_h2, total_derived;

                // if (p.UNPHASED)
                // {
                //     // unphased:
                //     //   num_1s = heterozygotes (one derived allele)
                //     //   num_2s = hom-derived genotypes (two derived alleles)
                //     derived_h1 = pass1_h1.num_1s[i_start] + 2 * pass1_h1.num_2s[i_start];
                //     derived_h2 = pass1_h2.num_1s[j_start] + 2 * pass1_h2.num_2s[j_start];
                // }
                // else
                // {
                //     // phased:
                //     //   num_1s already counts derived haplotypes
                //     derived_h1 = pass1_h1.num_1s[i_start];
                //     derived_h2 = pass1_h2.num_1s[j_start];
                // }

                // total_derived = derived_h1 + derived_h2;

                // // monomorphic if all alleles are ancestral or all are derived
                // if (total_derived == 0 || total_derived == total_alleles_joint)
                // {
                //     new_skip_h1.push_back(i_start);
                //     new_skip_h2.push_back(j_start);
                // }
            }
        }

        // Rebuild queues safely with old + new skips
        rebuild_skipqueue_with_new_skips(pass1_h1.skiplist, new_skip_h1);
        rebuild_skipqueue_with_new_skips(pass1_h2.skiplist, new_skip_h2);
    }

    else{
        // --------------------------------------------------------
        // Step 1: detect loci that become monomorphic when
        // combining hapData and hapData2
        // --------------------------------------------------------
        vector<size_t> new_skip_h1;
        vector<size_t> new_skip_h2;

        for (int i = 0; i < pass1_h1.nloci_before_filtering; i++) {

            int total_c1   = pass1_h1.num_1s[i] + pass1_h2.num_1s[i];
            int total_haps = pass1_h1.current_ngts + pass1_h2.current_ngts;

            // monomorphic if all alleles are ancestral or all derived
            if (total_c1 == 0 || total_c1 == total_haps) {
                new_skip_h1.push_back(i);
                new_skip_h2.push_back(i); // two copy for simple logic
                //cout<<"Locus " << i << " is monomorphic in combined data. Will skip from both files." << endl;
            }
        }
        rebuild_skipqueue_with_new_skips(pass1_h1.skiplist, new_skip_h1);
        rebuild_skipqueue_with_new_skips(pass1_h2.skiplist, new_skip_h2);
    }

    if(pass1_h2.nloci_before_filtering - pass1_h2.skiplist.size() == 0){
        HANDLE_ERROR("After filtering, no loci remain in reference VCF. Please check your input files.");
    }else if(pass1_h1.nloci_before_filtering - pass1_h1.skiplist.size() == 0){
        HANDLE_ERROR("After filtering, no loci remain in query VCF. Please check your input files.");
    }else{
        LOG("After filtering, " << pass1_h1.nloci_before_filtering - pass1_h1.skiplist.size() << " loci remain in query VCF.");
        LOG("After filtering, " << pass1_h2.nloci_before_filtering - pass1_h2.skiplist.size() << " loci remain in reference VCF.");
    }

    pass1_h1.skipcount = pass1_h1.skiplist.size();
    pass1_h2.skipcount = pass1_h2.skiplist.size();
    readHapDataVCF_pass2(filename, hapData, pass1_h1);
    readHapDataVCF_pass2(filename2, hapData2, pass1_h2);

}

void HapMap::readHapDataVCF(string filename, HapData& hapData)
{
    // RELATED FLAGS:
    // ARG_SKIP, ARG_KEEP, ARG_UNPHASED, ARG_LOW_MEM, ARG_MAF, ARG_MISSING

    initParamsInHap(hapData);

    igzstream fin;
    queue<size_t> skiplist;   // keep your original queue-based skip structure

    LOG("Opening " << filename << " to count number of haplotypes and loci...");

    fin.open(filename.c_str());
    if (fin.fail()) {
        HANDLE_ERROR("Failed to open " + filename + " for reading.");
    }

    constexpr int numMapCols = 9;   // standard fixed VCF columns before sample GT fields
    string line;

    size_t nloci_before_filtering = 0;   // number of variant records seen
    size_t skipcount = 0;                // number of records skipped
    int current_ngts = -1;               // number of genotype/sample columns

    // duplicate-position tracking
    int physpos = -1;
    int physpos_first_duplicated_id = -1;
    int skip_due_to_duplicate_pos = 0;

    // ------------------------------------------------------------
    // Small helpers for fast manual parsing
    // ------------------------------------------------------------

    // Return pointer to the end of current tab-delimited field
    auto find_field_end = [](const char* p, const char* end) -> const char* {
        while (p < end && *p != '\t') ++p;
        return p;
    };

    // Skip current field and land on first char of next field
    auto skip_to_next_field = [](const char* p, const char* end) -> const char* {
        while (p < end && *p != '\t') ++p;
        if (p < end) ++p;
        return p;
    };

    // Fast integer parse for POS column
    auto parse_int_fast = [](const char* s, const char* e) -> int {
        int x = 0;
        for (const char* p = s; p < e; ++p) {
            x = x * 10 + (*p - '0');
        }
        return x;
    };

    // Determine number of sample columns from #CHROM header line
    auto count_samples_from_header = [&](const string& header_line) -> int {
        int tabs = 0;
        for (char c : header_line) {
            if (c == '\t') ++tabs;
        }
        // total columns = tabs + 1, genotype columns = total - 9
        return (tabs + 1) - numMapCols;
    };

    // True only for single-base SNP alleles A/C/G/T
    auto is_snp_base = [](const char* s, const char* e) -> bool {
        if (e - s != 1) return false;
        const char c = *s;
        return c == 'A' || c == 'C' || c == 'G' || c == 'T';
    };

    // ============================================================
    // PASS 1:
    //   - find number of samples
    //   - count loci
    //   - decide which loci to skip
    //   - fill skiplist queue
    // ============================================================
    while (getline(fin, line))
    {
        if (line.empty()) continue;

        // Skip metadata/header lines. Parse #CHROM once to get sample count.
        if (line[0] == '#') {
            if (line.rfind("#CHROM", 0) == 0) {
                current_ngts = count_samples_from_header(line);
                if (current_ngts <= 0) {
                    HANDLE_ERROR("Failed to determine number of samples from VCF header in " + filename);
                }
            }
            continue;
        }

        if (current_ngts < 0) {
            HANDLE_ERROR("VCF header (#CHROM line) not found before variant data in " + filename);
        }

        nloci_before_filtering++;

        const char* pcur = line.c_str();
        const char* end  = pcur + line.size();

        // Parse first 9 fixed VCF columns once.
        const char* field_starts[numMapCols];
        const char* field_ends[numMapCols];

        for (int col = 0; col < numMapCols; ++col) {
            field_starts[col] = pcur;
            field_ends[col]   = find_field_end(pcur, end);
            pcur = field_ends[col];
            if (pcur < end) ++pcur;
        }

        // POS column
        const int new_physpos = parse_int_fast(field_starts[1], field_ends[1]);

        // Duplicate-position logic
        if (physpos == new_physpos) {
            skip_due_to_duplicate_pos += 1;
        } else {
            skip_due_to_duplicate_pos = 0;
            physpos_first_duplicated_id = static_cast<int>(nloci_before_filtering - 1);
        }
        physpos = new_physpos;

        bool skip_due_to_missing = false;
        bool skip_due_to_multiallelic = false;

        int number_of_1s = 0;
        int number_of_2s = 0;

        // If MULTI_ALLELIC mode is enabled, keep only simple SNPs.
        // This matches your earlier logic of skipping non-biallelic/non-SNP sites.
        if (p.MULTI_ALLELIC) {
            const char* ref_s = field_starts[3];
            const char* ref_e = field_ends[3];
            const char* alt_s = field_starts[4];
            const char* alt_e = field_ends[4];

            if (!is_snp_base(ref_s, ref_e) || !is_snp_base(alt_s, alt_e)) {
                skip_due_to_multiallelic = true;
            }
        }

        // Parse genotype/sample fields.
        // We only inspect the leading GT part, assuming fields look like:
        //   0|1
        //   0/1
        //   0|1:...
        //   1/1:...
        //
        // Like your original code, this does NOT fully parse multi-digit allele IDs.
        int seen_ngts = 0;

        if (p.UNPHASED)
        {
            // Unphased mode: treat genotype as diploid state per sample
            while (pcur < end && seen_ngts < current_ngts)
            {
                const char* g = pcur;

                // Need at least GT chars at positions 0,1,2
                if (end - g >= 3) {
                    const char allele1   = g[0];
                    const char separator = g[1];
                    const char allele2   = g[2];

                    // Missing / bad alleles
                    if (!p.MISSING_ALLOWED) {
                        if ((allele1 != '0' && allele1 != '1') ||
                            (allele2 != '0' && allele2 != '1'))
                        {
                            if (!p.MULTI_ALLELIC) {
                                HANDLE_ERROR("Alleles must be coded 0 or 1 only. Found alleles "
                                             << allele1 << separator << allele2);
                            } else {
                                skip_due_to_missing = true;
                            }
                        }
                    }

                    // If UNPHASED is on, '/' is allowed. If not, this still accepts '|'.
                    // No extra action needed here except for malformed alleles above.

                    // Branch-light genotype decode:
                    // '0' -> 0, '1' -> 1
                    const int a1 = (allele1 == '1');
                    const int a2 = (allele2 == '1');
                    const int sum = a1 + a2;

                    if (sum == 2) {
                        number_of_2s++;
                    } else if (sum == 1) {
                        number_of_1s++;
                    }
                } else {
                    if (!p.MULTI_ALLELIC) {
                        HANDLE_ERROR("Malformed genotype field at pos " << physpos);
                    } else {
                        skip_due_to_missing = true;
                    }
                }

                // Move to next sample field
                while (pcur < end && *pcur != '\t') ++pcur;
                if (pcur < end) ++pcur;
                seen_ngts++;
            }
        }
        else
        {
            // Phased mode: count haplotypes separately and reject '/' unless allowed by your old logic
            while (pcur < end && seen_ngts < current_ngts)
            {
                const char* g = pcur;

                if (end - g >= 3) {
                    const char allele1   = g[0];
                    const char separator = g[1];
                    const char allele2   = g[2];

                    if (!p.MISSING_ALLOWED) {
                        if ((allele1 != '0' && allele1 != '1') ||
                            (allele2 != '0' && allele2 != '1'))
                        {
                            if (!p.MULTI_ALLELIC) {
                                HANDLE_ERROR("Alleles must be coded 0 or 1 only. Found alleles "
                                             << allele1 << separator << allele2);
                            } else {
                                skip_due_to_missing = true;
                            }
                        }
                    }

                    if (separator != '|') {
                        if (!p.MULTI_ALLELIC) {
                            HANDLE_ERROR("Unphased entries detected (/ is used). Make sure you run with --unphased flag for correct results.");
                        } else {
                            skip_due_to_missing = true;
                        }
                    }

                    const int a1 = (allele1 == '1');
                    const int a2 = (allele2 == '1');

                    number_of_1s += a1 + a2;
                } else {
                    if (!p.MULTI_ALLELIC) {
                        HANDLE_ERROR("Malformed genotype field at pos " << physpos);
                    } else {
                        skip_due_to_missing = true;
                    }
                }

                while (pcur < end && *pcur != '\t') ++pcur;
                if (pcur < end) ++pcur;
                seen_ngts++;
            }
        }

        if (seen_ngts != current_ngts) {
            HANDLE_ERROR("line " + std::to_string(nloci_before_filtering) +
                         " of " + filename + " has " + std::to_string(seen_ngts) +
                         " genotype fields, but header indicates " + std::to_string(current_ngts) + ".");
        }

        const int nalleles_per_loc = current_ngts * 2;
        bool skip_due_to_maf = shouldSkipLocus(number_of_1s, number_of_2s, nalleles_per_loc);

        // Preserve your existing XP/XPNSL behavior
        if (p.CALC_XP || p.CALC_XPNSL) {
            skip_due_to_maf = false;
        }

        // If this is the second occurrence of a duplicate position,
        // also push the first occurrence to the skiplist.
        if (skip_due_to_duplicate_pos == 1) {
            if (skiplist.empty() || skiplist.back() != static_cast<size_t>(physpos_first_duplicated_id)) {
                if (physpos_first_duplicated_id != static_cast<int>(nloci_before_filtering - 2)) {
                    HANDLE_ERROR("Logic error in duplicate position handling.");
                }
                skiplist.push(nloci_before_filtering - 2);
                skipcount++;
            }
        }

        // Push current locus if any filter says skip
        if (skip_due_to_maf || skip_due_to_multiallelic || skip_due_to_missing || skip_due_to_duplicate_pos) {
            skiplist.push(nloci_before_filtering - 1);
            skipcount++;
        }
    }

    fin.clear();
    fin.close();

    // ============================================================
    // PASS 2:
    //   - allocate HapData exactly
    //   - reread file
    //   - skip loci in skiplist
    //   - load alleles into HapData
    // ============================================================
    fin.open(filename.c_str());
    if (fin.fail()) {
        HANDLE_ERROR("Failed to open " + filename + " for reading.");
    }

    if (p.SKIP && !p.CALC_XP && !p.CALC_XPNSL) {
        LOG(ARG_SKIP << " set. Removing all variants < " << p.MAF << ".\n");
    }

    const int nhaps = p.UNPHASED ? current_ngts : current_ngts * 2;

    LOG("Loading " << nhaps
                   << " haplotypes and " << (nloci_before_filtering - skipcount)
                   << " loci. Skipped " << skipcount << " loci \n");

    hapData.initHapData(nhaps, nloci_before_filtering - skipcount);

    // Preserve your old behavior
    hapData.skipQueue = skiplist;

    size_t locus = 0;
    size_t nloci_after_filtering = 0;

    while (getline(fin, line))
    {
        if (line.empty()) continue;
        if (line[0] == '#') continue;

        // Skip loci identified in pass 1
        if (!skiplist.empty() && skiplist.front() == locus) {
            skiplist.pop();
            locus++;
            continue;
        }

        const char* pcur = line.c_str();
        const char* end  = pcur + line.size();

        // Skip first 9 fixed VCF columns and land on first genotype field
        for (int col = 0; col < numMapCols; ++col) {
            pcur = skip_to_next_field(pcur, end);
        }

        if (p.UNPHASED)
        {
            // In unphased mode, one sample contributes one diploid state
            for (int field = 0; field < current_ngts; ++field)
            {
                // Branch-light decode of GT
                const int a1 = (pcur[0] == '1');
                const int a2 = (pcur[2] == '1');
                const int sum = a1 + a2;

                if (sum == 2) {
                    hapData.addAllele2(nloci_after_filtering, field);
                } else if (sum == 1) {
                    hapData.addAllele1(nloci_after_filtering, field);
                }

                while (pcur < end && *pcur != '\t') ++pcur;
                if (pcur < end) ++pcur;
            }
        }
        else
        {
            // In phased mode, each sample contributes two haplotypes
            for (int field = 0; field < current_ngts; ++field)
            {
                const int a1 = (pcur[0] == '1');
                const int a2 = (pcur[2] == '1');

                if (a1) hapData.addAllele1(nloci_after_filtering, 2 * field);
                if (a2) hapData.addAllele1(nloci_after_filtering, 2 * field + 1);

                while (pcur < end && *pcur != '\t') ++pcur;
                if (pcur < end) ++pcur;
            }
        }

        nloci_after_filtering++;
        locus++;
    }

    if (p.SKIP) {
        LOG("Removed " << skipcount << " low frequency variants from haplotype data.");
    }

    LOG("====");
    fin.close();
}

// Decide whether to skip a locus based on MAF and missing data
bool HapMap::shouldSkipLocus(int number_of_1s, int number_of_2s, int nalleles_per_loc)
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
    hapData.skipQueue = queue<size_t>();

    
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
            // readHapDataVCF(p.vcfFilename, *hapData); //query_vcf
            // readHapDataVCF(p.vcfFilename2, *hapData2); //ref_vcf
            readHapDataVCFXP(p.vcfFilename, p.vcfFilename2,  *hapData, *hapData2); // read both vcfs together and determine skiplists based on both
            // if (hapData->nloci != hapData2->nloci)
            // {
            //     HANDLE_ERROR("Haplotypes from " + p.vcfFilename + " and " + p.vcfFilename2 + " do not have the same number of loci.");
            // }
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