
#include "selscan-data.h"
#include <set>


void HapMap::initParamsInHap(HapData &hapData){
    hapData.MISSING_ALLOWED = p.MISSING_ALLOWED;
    hapData.unphased = p.UNPHASED;
}

void HapMap::readHapDataMSHap(string filename, HapData &hapData)
{
    initParamsInHap(hapData);
}


// START BITSET
//reads in haplotype data and also does basic checks on integrity of format
//returns a populated HaplotypeData structure if successful
//throws an exception otherwise
//impute hap IMPUTE HAP is transposed format where row represents loci,  column replesent individual
//so wc -l of impute hap is same as map.
void HapMap::readHapData(string filename, HapData& hapData)
{
    initParamsInHap(hapData);
    //PHASE 1: Read VCF File to get "nloci", "nhaps" and "skiplist"
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
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

        pair<int, int> fo = countFieldsAndOnes(line);
        current_nhaps = fo.first;
        //fo.secons is the number of 1s in the line (in locus)
        if( p.SKIP && (fo.second*1.0/current_nhaps < p.MAF || 1-(fo.second*1.0/current_nhaps) < p.MAF ) ) {
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

    //PHASE 2: Open VCF File To Load into Data Structure
    fin.open(filename.c_str());
    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    if(p.SKIP) { //prefilter all sites < MAF
        cerr << ARG_SKIP << " set. Removing all variants < " << p.MAF << ".\n";
        //(*flog)  << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
    }else{
        cerr << ARG_KEEP << " set. Not removing variants < " << p.MAF << ".\n";
        //(*flog) << ARG_KEEP << " set. Not removing variants < " << MAF << ".\n";
    }
    cerr << "Loading " << current_nhaps << " haplotypes and " << nloci_before_filter - skiplist.size()<< " loci...\n";

    if (p.UNPHASED){
        hapData.initHapData(current_nhaps/2, nloci_before_filter-skiplist.size());
    }else{
        hapData.initHapData(current_nhaps, nloci_before_filter-skiplist.size());
    }

    hapData.skipQueue = queue<int>();
    
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
        if (p.UNPHASED){
            line.erase(std::remove_if(line.begin(), line.end(), [](char c) {
                return std::isspace(static_cast<unsigned char>(c));
            }), line.end());

            if (current_nhaps % 2 != 0)
            {
                cerr << "ERROR:  Number of haplotypes must be even for unphased.\n";
                throw 0;
            }

            if(p.LOW_MEM){
                hapData.hapEntries[locus_after_filter].xorbitset->num_1s = 0;
                hapData.hapEntries[locus_after_filter].hapbitset->num_1s = 0;
            }
            
            for (int hap = 0; hap < current_nhaps; hap++){ 
                //xorbitset holds both 1 and 2
                if (hap % 2 == 0){
                    if(line[hap]=='1'  && line[hap+1]=='1'){
                        if(p.LOW_MEM){
                            hapData.hapEntries[locus_after_filter].xorbitset->set_bit(hap/2); // 2
                            hapData.hapEntries[locus_after_filter].xorbitset->num_1s += 1;
                        }
                    }else if ( (line[hap]=='1' && line[hap+1]=='0') || (line[hap]=='0' && line[hap+1]=='1') ){ // ==1
                        if(p.LOW_MEM){
                            hapData.hapEntries[locus_after_filter].hapbitset->set_bit(hap/2); // 1
                            hapData.hapEntries[locus_after_filter].hapbitset->num_1s += 1;
                        }
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
                    if(p.LOW_MEM){
                        hapData.hapEntries[locus_after_filter].hapbitset->set_bit(hap);
                        hapData.hapEntries[locus_after_filter].hapbitset->num_1s++;
                    }
                }
            }
        }

        locus_after_filter++;
        getline(fin, line);
    }
    fin.close();

    //PHASE 3: XOR and FLIP (flip disabled)
    hapData.xor_for_phased_and_unphased();
}


//tested-> unphased - lowmem, phased - lowmem
void HapMap::readHapDataVCF(string filename, HapData& hapData)
{
    initParamsInHap(hapData);
    // MULTI-chromosome support added
    set<string> chr_set;
    if(p.CHR_LIST!=""){
        std::stringstream ss(p.CHR_LIST);
        std::string item;
        while (std::getline(ss, item, ',')) { // Split the input string by commas
            if(item!="")
                chr_set.insert(item);
        }
    }

    if(p.MISSING_ALLOWED){
        cerr<<"Missing entries allowed: will impute missing entries if less than threshold"<<endl;
        readHapDataVCFMissing(filename, hapData);
        return;
    }
    //RELATED FLAGS: ARG_SKIP, ARG_KEEP, ARG_UNPHASED, ARG_LOW_MEM, ARG_MAF, ARG_MISSING
    
    
    igzstream fin;

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
        string chr;
        char allele1, allele2, separator;
        std::stringstream ss(line);
        int number_of_1s = 0;
        int number_of_2s = 0;

        for (int i = 0; i < numMapCols; i++) {
            ss >> junk;
            if(i==0){
                chr = junk;
                cout<<chr<<endl;
            }
        }
        for (int field = 0; field < current_nhaps; field++)
        {
            ss >> junk;
            allele1 = junk[0];
            separator = junk[1];
            allele2 = junk[2];

            if(!p.MISSING_ALLOWED){
                if ( (allele1 != '0' && allele1 != '1') || (allele2 != '0' && allele2 != '1') )
                {
                    cerr << "ERROR: Alleles must be coded 0/1 only.\n";
                    cerr << allele1 << " " << allele2 << endl;
                    throw 0;
                    exit(1);
                }
            }


            if(separator != '|' && !p.UNPHASED){
               cerr << "WARNING: Unphased entries detected. Make sure you run with --unphased flag for correct results.\n";
               throw 0;
            }

            if(p.UNPHASED){
                char allele = '0';
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

        int derived_allele_count = (p.UNPHASED? (number_of_1s + number_of_2s*2) : number_of_1s);

        if(!p.MISSING_ALLOWED){
            // --skip-low-freq filtering based on MAF
            bool skipreason1 = (p.SKIP && (derived_allele_count*1.0/(current_nhaps*2) < p.MAF || 1-(derived_allele_count*1.0/(current_nhaps*2)) < p.MAF ));
            bool skipreason2 = (chr_set.find(chr) != chr_set.end() && !chr_set.empty());
            if ( skipreason1 ||  skipreason2) {
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
        //(*flog)  << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
    }
   
    int nhaps = p.UNPHASED ? (current_nhaps ) : (current_nhaps ) * 2;
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

        for (int field = 0; field <  current_nhaps ; field++)
        {
            fin >> junk;
            allele1 = junk[0];
            separator = junk[1];

            allele2 = junk[2];
            if ( (allele1 != '0' && allele1 != '1') || (allele2 != '0' && allele2 != '1') )
            {
                cerr << "WARNING: Missing entry. ";
                cerr << allele1 << " and " << allele2 << endl;
                //throw 0;
            }

            if(p.UNPHASED){
                char allele = '0';
                if (allele1 == '1' && allele2 == '1'){
                    if(p.LOW_MEM){
                        hapData.hapEntries[nloci_after_filtering].xorbitset->set_bit(field);
                        hapData.hapEntries[nloci_after_filtering].xorbitset->num_1s++;
                    }
                }
                else if (allele1 == '1' || allele2 == '1'){
                    if(p.LOW_MEM){
                        hapData.hapEntries[nloci_after_filtering].hapbitset->set_bit(field);
                        hapData.hapEntries[nloci_after_filtering].hapbitset->num_1s++;
                    }
                }
            }else{ // phased
                if(allele1 == '1'){
                    if(p.LOW_MEM){
                        hapData.hapEntries[nloci_after_filtering].hapbitset->set_bit(2 * field);
                        hapData.hapEntries[nloci_after_filtering].hapbitset->num_1s++;
                    }
                }
                if(allele2 == '1'){
                    if(p.LOW_MEM){
                        hapData.hapEntries[nloci_after_filtering].hapbitset->set_bit(2 * field + 1);
                        hapData.hapEntries[nloci_after_filtering].hapbitset->num_1s++;
                    }
                }   
            }
        }
        nloci_after_filtering++;
    }

    // PHASE 3:  XOR
    hapData.xor_for_phased_and_unphased();
    

    

    if(p.SKIP){
        cerr << "Removed " << skipcount << " low frequency variants.\n";
        //(*flog) << "Removed " << skipcount << " low frequency variants.\n";
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
        //(*flog)  << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
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
                    if(p.LOW_MEM){
                        hapData.hapEntries[nloci_after_filtering].xorbitset->set_bit(field);
                        hapData.hapEntries[nloci_after_filtering].xorbitset->num_1s++;
                    }
                }
                else if ((allele1 == '1' && allele2 == '0') || (allele2 == '1' && allele1 == '0')){
                    if(p.LOW_MEM){
                        hapData.hapEntries[nloci_after_filtering].hapbitset->set_bit(field);
                        hapData.hapEntries[nloci_after_filtering].hapbitset->num_1s++;
                    }
                }else if(allele1 == '?' || allele2 == '?'){
                    if(p.LOW_MEM){
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
                                hapData.hapEntries[nloci_after_filtering].missbitset->set_bit(field);
                                hapData.hapEntries[nloci_after_filtering].missbitset->num_1s+= 1;
                            }
                        }
                    }
                }
            }else{ // phased
                if(allele1 == '1'){
                    if(p.LOW_MEM){
                        hapData.hapEntries[nloci_after_filtering].hapbitset->set_bit(2 * field);
                        hapData.hapEntries[nloci_after_filtering].hapbitset->num_1s++;
                    }
                }
                if(allele2 == '1'){
                    if(p.LOW_MEM){
                        hapData.hapEntries[nloci_after_filtering].hapbitset->set_bit(2 * field + 1);
                        hapData.hapEntries[nloci_after_filtering].hapbitset->num_1s++;
                    }
                }   
                if(allele1 == '?' || allele2 == '?'){
                    if(p.LOW_MEM){
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
                                hapData.hapEntries[nloci_after_filtering].xorbitset->set_bit(2 * field);
                                hapData.hapEntries[nloci_after_filtering].xorbitset->set_bit(2 * field + 1);
                                hapData.hapEntries[nloci_after_filtering].xorbitset->num_1s+= 2;
                            }
                        }
                    }
                }
            }
        }
        nloci_after_filtering++;
    }

    if(p.SKIP){
        cerr << "Removed " << skipcount << " low frequency variants.\n";
        //(*flog) << "Removed " << skipcount << " low frequency variants.\n";
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
        throw 0;
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
            throw 0;
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

        current_nhaps = p.UNPHASED? current_nhaps/2: current_nhaps; // each haplotype is represented by 2 columns
        for (int hap = 0; hap < current_nhaps; hap++)
        {
            if(p.UNPHASED){
                ss >> allele1;
                ss >> allele2;

                if ( (allele1 != '0' && allele1 != '1') || (allele2 != '0' && allele2 != '1') ){
                    cerr << "ERROR: Alleles must be coded 0/1 only. (Missing support not included for TPED)\n";
                    cerr << allele1 << " " << allele2 << endl;
                    throw 0;
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
                    throw 0;
                }
                if(allele=='1'){
                    number_of_1s++;
                }
            }
        }
        int derived_allele_count = (p.UNPHASED? (number_of_1s + number_of_2s*2) : number_of_1s);
        int total_allele_count = current_nhaps; // num rows
        current_nhaps = (p.UNPHASED? (current_nhaps/2) : (current_nhaps));
        if(!p.MISSING_ALLOWED){
            // --skip-low-freq filtering based on MAF
            bool skipreason1 = (p.SKIP && (derived_allele_count*1.0/(total_allele_count) < p.MAF || 1-(derived_allele_count*1.0/(total_allele_count)) < p.MAF ));
            bool skipreason2 = false;
            bool skipreason3 = false; //if(MISSING_ALLOWED){add here
            //bool skipreason2 = (chr_set.find(chr) != chr_set.end() && !chr_set.empty());
            if ( skipreason1 ||  skipreason2) {
                skipQueue_local.push(nloci_before_filter-1);
            }
        }
    }

    fin.clear(); // clear error flags
    //fin.seekg(fileStart);
    fin.close();

    int nloci_after_filter = nloci_before_filter-skipQueue_local.size();
    hapData.initHapData(current_nhaps, nloci_after_filter); // works for both phased and unphased
    cout<<"init by"<<hapData.nloci<<" "<<hapData.nhaps<<endl;
    // PHASE 2: Open again and load the data to HapData: now that we gathered some idea about the format of the file
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    cerr << "Loading " << current_nhaps << " haplotypes and " << nloci_after_filter << " loci...\n";
    if(p.SKIP){
        cerr << "Removed " << hapData.skipQueue.size() << " low frequency sites.\n";
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
                    hapData.hapEntries[locus].xorbitset->set_bit(hap); 
                    hapData.hapEntries[locus].xorbitset->num_1s++;
                }
                else if (allele1 == '1' || allele2 == '1'){ //allele = '1';
                    hapData.hapEntries[locus].hapbitset->set_bit(hap);
                    hapData.hapEntries[locus].hapbitset->num_1s++;
                }
            }else{
                char allele;
                fin >> allele;
                
                if(allele=='1'){
                    hapData.hapEntries[locus].hapbitset->set_bit(hap);
                    hapData.hapEntries[locus].hapbitset->num_1s++;
                }
            }
        }
    }

    fin.close();

    //print(5);
}

