
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
    //         exit(2);
    //     }
    // }

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
    initParamsInHap(hapData);

    //PHASE 1: Read input file to get "nloci", "nhaps" and "skiplist"
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        exit(EXIT_FAILURE);
    }

    string line, prev_line;
    int current_nhaps = 0;
    int num_rows = 0;

    queue<int> skiplist;

    getline(fin, prev_line); // first line
    int nloci_before_filter =  countFields(prev_line);
    vector<int> number_of_1s_per_locus(nloci_before_filter, 0);
    vector<int> number_of_2s_per_locus(nloci_before_filter, 0);

    //Counts number of haps (rows) and number of loci (cols)
    //if any lines differ, send an error message and throw an exception

    int hap = 1;
    while (getline(fin, line))
    {
        if(prev_line.length() != line.length()){
            cerr << "ERROR: line " << hap << " of " << filename << " has " <<  line.length()/2
                 << " lines, but the previous line has " << prev_line.length()/2 << ".\n";
            exit(EXIT_FAILURE);
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
                if(line[loc]=='1'){
                    number_of_1s_per_locus[loc]++;
                }
            }
        }
            
        prev_line = line; // store the previous line
        hap++;
    }
    
    num_rows = hap;
    current_nhaps = (p.UNPHASED) ? hap/2 : hap;

    for(int loc=0; loc<nloci_before_filter; loc++){
        if(skipThisLocus(number_of_1s_per_locus[loc], number_of_2s_per_locus[loc], num_rows)){
            skiplist.push(loc);
        }
    }
    
    fin.clear(); // clear error flags
    fin.close();

    number_of_1s_per_locus.clear();
    number_of_2s_per_locus.clear();

    if(p.SKIP) { //prefilter all sites < MAF
        cerr << ARG_SKIP << " set. Removing all variants < " << p.MAF << ".\n";
        //(*flog)  << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
    }else{
        cerr << ARG_KEEP << " set. Not removing variants < " << p.MAF << ".\n";
        //(*flog) << ARG_KEEP << " set. Not removing variants < " << MAF << ".\n";
    }
    
    //current_nhaps = (p.UNPHASED) ? current_nhaps/2 : current_nhaps;
    
    cerr << "Loading " << current_nhaps << " haplotypes and " << nloci_before_filter - skiplist.size()<< " loci...\n";
    hapData.initHapData(current_nhaps, nloci_before_filter-skiplist.size());
    hapData.skipQueue = queue<int>();

    //PHASE 2: Open input File 2nd time To Load into Data Structure
    fin.open(filename.c_str());
    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        exit(EXIT_FAILURE);
    }
    //getline(fin, line);
    stringstream ss;
    char allele1;
    int locus_after_filter = 0;
    for (int hap = 0; hap < num_rows; hap++)
    {
        // if(!skiplist.empty()){
        //     if(skiplist.front() == locus){
        //         skiplist.pop();
        //         hapData.skipQueue.push(locus);
        //         getline(fin, line);
        //         continue;
        //     }
        // }
        getline(fin, line);
        ss.str(line);
        line.erase(std::remove_if(line.begin(), line.end(), [](char c) {
                return std::isspace(static_cast<unsigned char>(c));
            }), line.end());


        for(int loc=0; loc<nloci_before_filter; loc++){
            if(p.UNPHASED){
                if (hap % 2 == 1){
                    if(line[loc]=='1'  && prev_line[loc]=='1'){
                        hapData.addAllele2(loc, (hap-1)/2);
                    }else if ( (line[loc]=='1' && prev_line[loc]=='0') || (line[loc]=='0' && prev_line[loc]=='1') ){ // ==1
                        hapData.addAllele1(loc, (hap-1)/2);
                    }
                }
            }else{
                if(line[loc]=='1'){
                    hapData.addAllele1(loc, hap);
                }else{
                    if(line[loc]!='0'){
                        cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
                        exit(EXIT_FAILURE);
                    }

                }
            }
        }
        getline(fin, line);
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
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        exit(EXIT_FAILURE);
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
        current_nhaps = countFields(line);
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

        if(skipThisLocus(number_of_1s, number_of_2s, current_nhaps)){
            skiplist.push(nloci_before_filter-1);
        }
    
        if (previous_nhaps < 0)
        {
            previous_nhaps = current_nhaps;
            continue;
        }
        else if (previous_nhaps != current_nhaps)
        {
            cerr << "ERROR: line " << nloci_before_filter << " of " << filename << " has " << current_nhaps
                 << ", but the previous line has " << previous_nhaps << ".\n";
            throw 0;
        }
        previous_nhaps = current_nhaps;
    }
    fin.clear(); // clear error flags
    fin.close();

    if(p.SKIP) { //prefilter all sites < MAF
        cerr << ARG_SKIP << " set. Removing all variants < " << p.MAF << ".\n";
        //(*flog)  << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
    }else{
        cerr << ARG_KEEP << " set. Not removing variants < " << p.MAF << ".\n";
        //(*flog) << ARG_KEEP << " set. Not removing variants < " << MAF << ".\n";
    }
    
    current_nhaps = (p.UNPHASED) ? current_nhaps/2 : current_nhaps;
    cerr << "Loading " << current_nhaps << " haplotypes and " << nloci_before_filter - skiplist.size()<< " loci...\n";
    hapData.initHapData(current_nhaps, nloci_before_filter-skiplist.size());
    hapData.skipQueue = queue<int>();

    //PHASE 2: Open input File 2nd time To Load into Data Structure
    fin.open(filename.c_str());
    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        exit(EXIT_FAILURE);
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

        if (p.UNPHASED){
            if (current_nhaps % 2 != 0)
            {
                cerr << "ERROR:  Number of haplotypes must be even for unphased.\n";
                throw 0;
            }
            
            for (int hap = 0; hap < current_nhaps; hap++){ 
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
                    cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
                    throw 0;
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
    initParamsInHap(hapData);
    if(p.MISSING_ALLOWED){
        cerr<<"WARNING: Missing entries allowed: will impute missing entries if less than threshold."<<endl;
        *flog<<"WARNING: Missing entries allowed: will impute missing entries if less than threshold."<<endl;

        readHapDataVCFMissing(filename, hapData);
        return;
    }
    //RELATED FLAGS: ARG_SKIP, ARG_KEEP, ARG_UNPHASED, ARG_LOW_MEM, ARG_MAF, ARG_MISSING
    
    igzstream fin;

    queue<int> skiplist;

    // PHASE 1: Counting so that inititalization is smooth
    cerr << "Opening " << filename << "...\n";
    *flog << "Opening " << filename << "...\n";

    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        *flog << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
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
                    cerr << "ERROR: Alleles must be coded 0 or 1 only. Found alleles "<< allele1 << separator << allele2 << endl;
                    *flog << "ERROR: Alleles must be coded 0 or 1 only. Found alleles "<< allele1 << separator << allele2 << endl;
                    throw 0;
                    exit(EXIT_FAILURE);
                }
            }

            if(separator != '|' && !p.UNPHASED){
               cerr << "ERROR: Unphased entries detected (| is used). Make sure you run with --unphased flag for correct results.\n";
               throw 0;
               exit(EXIT_FAILURE);
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
        bool skipreason1 = skipThisLocus(number_of_1s, number_of_2s, nalleles_per_loc);
        bool skipreason2 = false;
        
        // if(p.MULTI_CHR){ 
        //     if(chr_set.empty()){
        //         //cerr<<"WARNING: No chromosome list provided. Running analysis on all chromosomes.\n";
        //     }else{
        //         skipreason2 = (chr_set.find(chr) != chr_set.end());
        //     }
        // }

        if(!p.CALC_XP &&  !p.CALC_XPNSL){
            if ( skipreason1 ||  skipreason2) {
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
            cerr << "ERROR: line " << nloci_before_filtering << " of " << filename << " has " << current_ngts
                 << " fields, but the previous line has " << prev_ngts << " fields.\n";
            throw 0;
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
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    if(p.SKIP && !p.CALC_XP &&  !p.CALC_XPNSL){
        cerr << ARG_SKIP << " set. Removing all variants < " << p.MAF << ".\n";
        //(*flog)  << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
    }
   
    int nhaps = p.UNPHASED ? (current_ngts ) : (current_ngts ) * 2;
    cerr << "Loading " << nhaps << " haplotypes and " << nloci_before_filtering-skipcount << " loci...\n";
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
            // if ( (allele1 != '0' && allele1 != '1') || (allele2 != '0' && allele2 != '1') )
            // {
            //     cerr << "WARNING: Missing entry. ";
            //     cerr << allele1 << " and " << allele2 << endl;
            //     //throw 0;
            // }

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
        cerr << "Removed " << skipcount << " low frequency variants.\n";
        (*flog) << "Removed " << skipcount << " low frequency variants.\n";
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
    //  exit(1);
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
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
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
            

            //if(separator != '|'){
            //    cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
            //    throw 0;
            //}

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
                    cerr << "WARNING: Alleles not coded 0/1, treating this as missing site.\n";
                    // cerr << allele1 << " " << allele2 << endl;
                    // throw 0;

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
                cout<<"Skipping cause Derived "<<derived_allele_count<<" Missing "<<missing_count<<" Total "<<current_nhaps*2<<endl;
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
            cerr << "ERROR: line " << nloci_before_filtering << " of " << filename << " has " << current_nhaps
                 << " fields, but the previous line has " << previous_nhaps << " fields.\n";
            throw 0;
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
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    if(p.SKIP){
        cerr << ARG_SKIP << " set. Removing all variants < " << p.MAF << ".\n";
        (*flog)  << ARG_SKIP << " set. Removing all variants < " << p.MAF << ".\n";
    }
   
    int nhaps = p.UNPHASED ? (current_nhaps ) : (current_nhaps ) * 2;
    cerr << "Loading " << nhaps << " haplotypes and " << nloci_before_filtering-skipcount << " loci... skipped "<<skipcount << "loci \n";
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
                cerr << "ERROR: Alleles must be coded 0/1 only.\n";
                cerr << allele1 << " " << allele2 << endl;

                allele1 = '?';
                allele2 = '?';
                //throw 0;
            }

            //if(separator != '|'){
            //    cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
            //    throw 0;
            //}

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
        cerr << "Removed " << skipcount << " low frequency variants.\n";
        (*flog) << "Removed " << skipcount << " low frequency variants.\n";
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
    //  exit(1);
}



bool HapMap::skipThisLocus(int number_of_1s, int number_of_2s, int nalleles_per_loc, int missing_count){ 
    //nalleles same for phased and unphased, tot n of columns

    int derived_allele_count = (p.UNPHASED? (number_of_1s + number_of_2s*2) : number_of_1s);
    
    return (derived_allele_count)*1.0/(nalleles_per_loc) == 0 ;
    return 1-(derived_allele_count*1.0/(nalleles_per_loc)) == 0;

    bool skipreason1 = false; // --skip-low-freq filtering based on MAF
    bool skipreason2 = false; // missing

    double ALLOWED_MISSING_PERC = 1;
    
    if(p.MISSING_ALLOWED){
        int ancestral_allele_count = nalleles_per_loc - derived_allele_count;
        bool derived_lower = (derived_allele_count+missing_count)*1.0/(nalleles_per_loc) < this->MIN_MAF_CUTOFF ;
        bool ancestral_lower = (ancestral_allele_count+missing_count)*1.0/(nalleles_per_loc) < this->MIN_MAF_CUTOFF ;
        if(p.SKIP){
                skipreason2 = (derived_lower || ancestral_lower);
        }

        bool missing_cutoff_crossed = 1.0*missing_count/nalleles_per_loc > ALLOWED_MISSING_PERC;
        skipreason2 = (skipreason2 || missing_cutoff_crossed);
        // if (missing_cutoff_crossed) {
        //     //cout<<"Skipping cause Derived "<<derived_allele_count<<" Missing "<<missing_count<<" Total "<<current_nhaps*2<<endl;
        //     skipreason2 = true;
        // }
    }else{
        skipreason1 = (p.SKIP && (derived_allele_count*1.0/(nalleles_per_loc) < this->MIN_MAF_CUTOFF || 1-(derived_allele_count*1.0/(nalleles_per_loc)) < this->MIN_MAF_CUTOFF ));

    }
    return (skipreason1 ||  skipreason2 );
}
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
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());
    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        exit(2);
    }

    int numMapCols = 4;
    //int fileStart = fin.tellg();
    string line;
    int nloci_before_filter = 0; // number of loci before filtering (filtering happens when SKIP is enabled)
    int previous_nhaps = -1;
    int current_nhaps = 0;
    
    queue<int> skipQueue_local;

    //Counts number of haps (cols) and number of loci (rows)
    //if any lines differ, send an error message and throw an exception
    while (getline(fin, line))
    {
        //getline(fin,line);
        //if(fin.eof()) break;
        nloci_before_filter++; // each line contributes to a locus
        current_nhaps = countFields(line) - numMapCols; // here, current_nhaps (count of fields) is the number of haps 
        
        if (previous_nhaps < 0)
        {
            previous_nhaps = current_nhaps;
            continue;
        }

        else if (previous_nhaps != current_nhaps)
        {
            cerr << "ERROR: line " << nloci_before_filter << " of " << filename << " has " << current_nhaps
                 << " fields, but the previous line has " << previous_nhaps << " fields.\n";
            exit(2);
        }
        previous_nhaps = current_nhaps;

        int number_of_1s = 0;
        int number_of_2s = 0;

        stringstream ss;
        ss.str(line);
          
        char allele1, allele2;
        string junk;
        
        for(int i = 0; i < numMapCols; i++){
            ss >> junk;
        }

        for (int hap = 0; hap < current_nhaps; hap++)
        {
            if(p.UNPHASED){
                ss >> allele1;
                ss >> allele2;

                if ( (allele1 != '0' && allele1 != '1') || (allele2 != '0' && allele2 != '1') ){
                    cerr << "ERROR: Alleles must be coded 0/1 only. (Missing support not included for TPED)\n";
                    cerr << allele1 << " " << allele2 << endl;
                    exit(2);
                }
                
                if (allele1 == '1' && allele2 == '1'){
                    number_of_2s++; //allele = '2';
                }
                else if (allele1 == '1' || allele2 == '1'){
                    number_of_1s++; //allele = '1';
                }
            }else{ //phased
                char allele;
                ss >> allele;
                if (allele!= '0' && allele!= '1')
                {
                    cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
                    cerr << "Allele "<<allele << endl;
                    exit(2);
                }
                if(allele=='1'){
                    number_of_1s++;
                }
            }
        }
        if(skipThisLocus(number_of_1s, number_of_2s, current_nhaps)){
            skipQueue_local.push(nloci_before_filter-1);
        }
    }

    fin.clear(); // clear error flags
    //fin.seekg(fileStart);
    fin.close();

    int nloci_after_filter = nloci_before_filter-skipQueue_local.size();
    current_nhaps = p.UNPHASED? current_nhaps/2: current_nhaps; // each haplotype is represented by 2 columns
    hapData.initHapData(current_nhaps, nloci_after_filter); // works for both phased and unphased
    
    
    // PHASE 2: Open again and load the data to HapData: now that we gathered some idea about the format of the file
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        exit(2);
    }

    cerr << "Loading " << current_nhaps << " haplotypes and " << nloci_after_filter << " loci...\n";
    if(p.SKIP){
        cerr << "Removed " << hapData.skipQueue.size() << " low frequency variants.\n";
    }

    string junk;
    for (int locus = 0; locus < nloci_before_filter; locus++)
    {
        if(!skipQueue_local.empty()){
            if(skipQueue_local.front() == locus){
                skipQueue_local.pop();    
                hapData.skipQueue.push(locus); // make a copy to skipQueue
                getline(fin, junk);
                continue;
            }
        }
        for (int i = 0; i < numMapCols; i++)
        {
            fin >> junk; // skip first 4 columns
        }
        for (int hap = 0; hap < hapData.nhaps; hap++)
        {
            if(p.UNPHASED){
                char allele1, allele2;
                fin >> allele1;
                fin >> allele2;
                
                if (allele1 == '1' && allele2 == '1'){ //allele = '2';
                    hapData.addAllele2(locus, hap);
                }
                else if (allele1 == '1' || allele2 == '1'){ //allele = '1';
                    hapData.addAllele1(locus, hap);
                }
            }else{
                char allele;
                fin >> allele;
                
                if(allele=='1'){
                    hapData.addAllele1(locus, hap);
                }
            }
        }
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
                std::cerr << "ERROR: Haplotypes from " << p.tpedFilename << " and " << p.tpedFilename2 << " do not have the same number of loci.\n";
                exit(2);
            }
        }
        mapData->readMapDataTPED(p.tpedFilename, hapData->nloci, hapData->nhaps, p.USE_PMAP, hapData->skipQueue);
    }
    else if (p.VCF) {
        if (p.CALC_XP || p.CALC_XPNSL)
        {
            readHapDataVCF(p.vcfFilename, *hapData);
            readHapDataVCF(p.vcfFilename2, *hapData2);
            if (hapData->nloci != hapData2->nloci)
            {
                std::cerr << "ERROR: Haplotypes from " << p.vcfFilename << " and " << p.vcfFilename2 << " do not have the same number of loci.\n";
                exit(2);
            }
        }else{
            readHapDataVCF(p.vcfFilename, *hapData);
        }
        if(!p.CALC_NSL && !p.CALC_XPNSL && !p.USE_PMAP) {
            mapData->readMapData(p.mapFilename, hapData->nloci, p.USE_PMAP, hapData->skipQueue);
        }
        else{//Load physical positions
            mapData->readMapDataVCF(p.vcfFilename, hapData->nloci, hapData->skipQueue);
        }
    }else if(p.THAP){
        readHapDataTHAP(p.hapFilename, *hapData);
        if (p.CALC_XP || p.CALC_XPNSL)
        {
            readHapDataTHAP(p.hapFilename2, *hapData2);
            if (hapData->nloci != hapData2->nloci)
            {
                std::cerr << "ERROR: Haplotypes from " << p.hapFilename << " and " << p.hapFilename2 << " do not have the same number of loci.\n";
                exit(2);
            }
        }
        mapData->readMapData(p.mapFilename, hapData->nloci, p.USE_PMAP, hapData->skipQueue);
    }else
    {
        readHapDataMSHap(p.thapFilename, *hapData);
        if (p.CALC_XP || p.CALC_XPNSL)
        {
            readHapDataMSHap(p.thapFilename2, *hapData2);
            if (hapData->nloci != hapData2->nloci)
            {
                std::cerr << "ERROR: Haplotypes from " << p.thapFilename << " and " << p.thapFilename2 << " do not have the same number of loci.\n";
                exit(2);
            }
        }
        mapData->readMapData(p.mapFilename, hapData->nloci, p.USE_PMAP, hapData->skipQueue);
    }


    // PHASE 2:  XOR
    if(p.benchmark_flag=="XOR")
        hapData->xor_for_phased_and_unphased();


    auto end_reading = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> read_duration =  end_reading - start_reading;
    std::cerr<<("LOG: Input file loaded in "+to_string(read_duration.count())+" s.")<<endl;
    (*flog)<<("LOG: Input file loaded in "+to_string(read_duration.count())+" s.\n")<<endl;;
    
    // DEBUG::: mapData.print();

    // if(mapData->chr_list.size() > 1 && p.MULTI_CHR == false){
    //     cerr<<"ERROR: Input file contains multiple chromosomes, but no chromosome list is provided. Use --multi-chr flag to specify chromosomes.\n";
    //     exit(2);
    // }

    // Check if map is in order
    for (int i = 1; i < mapData->nloci; i++) {
        if ( mapData->mapEntries[i].physicalPos <  mapData->mapEntries[i-1].physicalPos ) {
            std::cerr << "ERROR: Variant physical position must be monotonically increasing.\n";
            std::cerr << "\t" << mapData->mapEntries[i].locusName << " " << mapData->mapEntries[i].physicalPos << " appears after";
            std::cerr << "\t" <<  mapData->mapEntries[i-1].locusName << " " << mapData->mapEntries[i-1].physicalPos << "\n";
            exit(2);
        }
        if ( !p.CALC_NSL && mapData->mapEntries[i].geneticPos  < mapData->mapEntries[i-1].geneticPos  ) {
            std::cerr << "ERROR: Variant genetic position must be monotonically increasing.\n";
            std::cerr << "\t" << mapData->mapEntries[i].locusName << " " << mapData->mapEntries[i].geneticPos << " appears after";
            std::cerr << "\t" << mapData->mapEntries[i-1].locusName << " " << mapData->mapEntries[i-1].geneticPos << "\n";
            exit(2);
        }
    }
    
}