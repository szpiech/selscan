#include "hapdata.h"
#include "../gzstream.h"
#include "../selscan-cli.h"
#include <algorithm>
#include <sstream>
#include <memory>
#include <limits>

HapData::~HapData(){
    releaseHapData();
}


/** Sets up structure according to nhaps and nloci
 * 
*/
void HapData::initHapData(int nhaps, int nloci)
{
    if (nhaps > std::numeric_limits<int>::max() || nloci > std::numeric_limits<int>::max())
    {
        cerr << "ERROR: number of haplotypes (" << nhaps << ") and number of loci (" << nloci << ") exceed maximum value.\n";
        throw 0;
    }

    if (nhaps < 1 || nloci < 1)
    {
        cerr << "ERROR: number of haplotypes (" << nhaps << ") and number of loci (" << nloci << ") must be positive.\n";
        throw 0;
    }
    cout << "Final number of haplotypes " << nhaps << ", number of loci " << nloci << ".\n";

    this->hapEntries = new struct HapEntry[nloci];
    this->nhaps = nhaps;
    this->nloci = nloci;

    if(LOW_MEM){
        for (int j = 0; j < nloci; j++){
            hapEntries[j].hapbitset = new MyBitset(nhaps);
            hapEntries[j].xorbitset = new MyBitset(nhaps);
        }
        INIT_SUCCESS = true;
    }else{
        INIT_SUCCESS = true;
    }
}

void HapData::releaseHapData()
{
    if(LOW_MEM){
        if (hapEntries == NULL) return;

        //we have a MyBitset for every locus
        for (int j = 0; j < nloci; j++){
            delete hapEntries[j].hapbitset ; //MyBitset destructor called
            delete hapEntries[j].xorbitset; //MyBitset destructor called
        }
        delete [] hapEntries;
        hapEntries = NULL;
        this->nhaps = -9;
        this->nloci = -9;
        return;
    }else{
        if (hapEntries == NULL) return;
        for (int i = 0; i < this->nloci; i++)
        {
            //delete [] hapEntries[i];
            //delete hapEntries[i];
            hapEntries[i].xors.clear();
            hapEntries[i].positions.clear();    
        }

        delete [] hapEntries;

        hapEntries = NULL;
        this->nhaps = -9;
        this->nloci = -9;
        //data = NULL;
        return;
    }
}

void HapData::xor_for_phased_and_unphased(){
    if(unphased){
        for(int nloci_after_filtering = 0; nloci_after_filtering < this->nloci; nloci_after_filtering++){
            if( !LOW_MEM ){
                getThreeUnphasedGroups(hapEntries[nloci_after_filtering].positions, hapEntries[nloci_after_filtering-1].positions,hapEntries[nloci_after_filtering].positions2, hapEntries[nloci_after_filtering-1].positions2, hapEntries[nloci_after_filtering].g);
            }
        }
    }else{
        if(benchmark_flag != "XOR"){
            return;
        } 
        for(int nloci_after_filtering = 0; nloci_after_filtering < this->nloci; nloci_after_filtering++){
            if(nloci_after_filtering==0){
                // if(nloci<=1){
                //     throw "ERROR: Dataset has only 1 locus, XOR out of bound.";    
                // }

                if(LOW_MEM){
                    // CHANGEXOR
                    MyBitset* b1 =(hapEntries[0].hapbitset);
                    MyBitset* b2 = (hapEntries[1].hapbitset);
                    int sum = 0;
                    for (int k = 0; k < b1->nwords; k++) {
                        hapEntries[0].xorbitset->bits[k] = b1->bits[k] ^ b2->bits[k];
                        sum += __builtin_popcountll(hapEntries[0].xorbitset->bits[k]);
                    }
                    hapEntries[0].xorbitset->num_1s = sum;
                }else{
                    //CHANGEXOR
                    std::set_symmetric_difference(hapEntries[0].positions.begin(), hapEntries[0].positions.end(),hapEntries[1].positions.begin(), hapEntries[1].positions.end(), std::back_inserter(hapEntries[0].xors));
                }
            }else{
                if(LOW_MEM){
                    MyBitset* b1 =(hapEntries[nloci_after_filtering].hapbitset);
                    MyBitset* b2 = (hapEntries[nloci_after_filtering-1].hapbitset);

                    int sum = 0;
                    for (int k = 0; k < b1->nwords; k++) {
                        hapEntries[nloci_after_filtering].xorbitset->bits[k] = b1->bits[k] ^ b2->bits[k];
                        sum += __builtin_popcountll(hapEntries[nloci_after_filtering].xorbitset->bits[k]);
                    }
                    hapEntries[nloci_after_filtering].xorbitset->num_1s = sum;
                }else{
                    vector<int>& curr_xor = hapEntries[nloci_after_filtering].xors;
                    vector<int>& curr_positions = hapEntries[nloci_after_filtering].positions;
                    vector<int>& prev_positions = hapEntries[nloci_after_filtering-1].positions;
                    std::set_symmetric_difference(curr_positions.begin(), curr_positions.end(),prev_positions.begin(), prev_positions.end(),
                                    std::back_inserter(curr_xor));
                }                    
            }
        }

    }    
}

// START BITSET
//reads in haplotype data and also does basic checks on integrity of format
//returns a populated HaplotypeData structure if successful
//throws an exception otherwise
//impute hap IMPUTE HAP is transposed format (thap) where row represents loci,  column replesent individual
//so wc -l of impute hap is same as map.
void HapData::readHapData(string filename)
{
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

    vector<int> num_1s_per_loci;

    //Counts number of haps (rows) and number of loci (cols)
    //if any lines differ, send an error message and throw an exception
    while (getline(fin, line))
    {
        nloci_before_filter++;

        pair<int, int> fo = countFieldsAndOnes(line);
        current_nhaps = fo.first;
        num_1s_per_loci.push_back(fo.second);
        if( SKIP && (fo.second*1.0/current_nhaps < MAF || 1-(fo.second*1.0/current_nhaps) < MAF ) ) {
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

    if(SKIP) { //prefilter all sites < MAF
        cerr << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
        (*flog)  << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
    }else{
        cerr << ARG_KEEP << " set. Not removing variants < " << MAF << ".\n";
        (*flog) << ARG_KEEP << " set. Not removing variants < " << MAF << ".\n";
    }
    cerr << "Loading " << current_nhaps << " haplotypes and " << nloci_before_filter - skiplist.size()<< " loci...\n";

    if (unphased){
        initHapData(current_nhaps/2, nloci_before_filter-skiplist.size());
    }else{
        initHapData(current_nhaps, nloci_before_filter-skiplist.size());
    }

    this->skipQueue = queue<int>();
    
    getline(fin, line);
    stringstream ss;
    char allele1;
    int locus_after_filter = 0;
    for (int locus = 0; locus < nloci_before_filter; locus++)
    {
        if(!skiplist.empty()){
            if(skiplist.front() == locus){
                skiplist.pop();
                skipQueue.push(locus);
                getline(fin, line);
                continue;
            }
        }

        ss.str(line);
        if (unphased){
            line.erase(std::remove_if(line.begin(), line.end(), [](char c) {
                return std::isspace(static_cast<unsigned char>(c));
            }), line.end());

            if (current_nhaps % 2 != 0)
            {
                cerr << "ERROR:  Number of haplotypes must be even for unphased.\n";
                throw 0;
            }

            if(LOW_MEM){
                this->hapEntries[locus_after_filter].xorbitset->num_1s = 0;
                this->hapEntries[locus_after_filter].hapbitset->num_1s = 0;
            }
            
            for (int hap = 0; hap < current_nhaps; hap++){ 
                //xorbitset holds both 1 and 2
                if (hap % 2 == 0){
                    if(line[hap]=='1'  && line[hap+1]=='1'){
                        if(LOW_MEM){
                            hapEntries[locus_after_filter].xorbitset->set_bit(hap/2); // 2
                            hapEntries[locus_after_filter].xorbitset->num_1s += 1;
                        }else{
                            hapEntries[locus_after_filter].positions2.push_back(hap/2);
                        }
                    }else if ( (line[hap]=='1' && line[hap+1]=='0') || (line[hap]=='0' && line[hap+1]=='1') ){ // ==1
                        if(LOW_MEM){
                            this->hapEntries[locus_after_filter].hapbitset->set_bit(hap/2); // 1
                            this->hapEntries[locus_after_filter].hapbitset->num_1s += 1;
                        }else{
                            hapEntries[locus_after_filter].positions.push_back(hap/2);
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
                    if(LOW_MEM){
                        this->hapEntries[locus_after_filter].hapbitset->set_bit(hap);
                        this->hapEntries[locus_after_filter].hapbitset->num_1s++;
                    }else{
                        this->hapEntries[locus_after_filter].positions.push_back(hap);
                    }
                }
            }
        }

        locus_after_filter++;
        getline(fin, line);
    }
    fin.close();


    //PHASE 3: XOR
    xor_for_phased_and_unphased();
    
    //PHASE 4: FLIP
    /*
    for (int locus_after_filter = 0; locus_after_filter < this->nloci; locus_after_filter++){
        if(hapEntries[locus_after_filter].hapbitset->num_1s > nhaps/2){
            hapEntries[locus_after_filter].flipped = true;
            MyBitset* b1;
            b1 = hapEntries[locus_after_filter].hapbitset;

            //#pragma omp simd
            for(int k = 0; k<b1->nwords; k++){
                b1->bits[k] = ~(b1->bits[k]);   // negate all bits
            }

            //#pragma omp simd
            for(int i = b1->nbits; i<b1->nwords*b1->WORDSZ; i++){
                b1->clear_bit(i);       // clear the trailing bits
            }
        }else{
            hapEntries[locus_after_filter].flipped = false;
        }
    }
    */


    // // // DEBUG
    // for(int locus_after_filter = 0; locus_after_filter < 1; locus_after_filter++){
    //     // cout<<locus_after_filter<<"::: ";
    //     // hapEntries[locus_after_filter].hapbitset->print_pos();

    //     cout<<locus_after_filter<<"::: ";
    //     hapEntries[locus_after_filter].xorbitset->print_pos();

    // }
    // exit(1);
}


//tested-> unphased - lowmem, phased - lowmem
void HapData::readHapDataVCF(string filename)
{
    if(MISSING){
        readHapDataVCFMissing(filename);
        return;
    }
    //RELATED FLAGS: ARG_SKIP, ARG_KEEP, ARG_UNPHASED, ARG_LOW_MEM, ARG_MAF, ARG_MISSING
    


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
        for (int field = 0; field < current_nhaps; field++)
        {
            ss >> junk;
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

        int derived_allele_count = (unphased? (number_of_1s + number_of_2s*2) : number_of_1s);

        // --skip-low-freq filtering based on MAF
        if ( SKIP && (derived_allele_count*1.0/(current_nhaps*2) < MAF || 1-(derived_allele_count*1.0/(current_nhaps*2)) < MAF ) ) {
            skiplist.push(nloci_before_filtering-1);
            skipcount++;
        } else {
            number_of_1s_per_loci.push_back(number_of_1s);
            if(unphased){
                number_of_2s_per_loci.push_back(number_of_2s);
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

    if(SKIP){
        cerr << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
        (*flog)  << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
    }
   
    int nhaps = unphased ? (current_nhaps ) : (current_nhaps ) * 2;
    cerr << "Loading " << nhaps << " haplotypes and " << nloci_before_filtering-skipcount << " loci...\n";
    initHapData(nhaps, nloci_before_filtering-skipcount);

    string junk;
    char allele1, allele2, separator;
    bool skipLine = false; // to skip metadata lines

    skipQueue = skiplist; 
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
        
        if(unphased){
            if(benchmark_flag == "XOR"){
                hapEntries[nloci_after_filtering].xors.reserve(number_of_1s_per_loci[nloci_after_filtering]+number_of_2s_per_loci[nloci_after_filtering]);
            }

            hapEntries[nloci_after_filtering].positions.reserve(number_of_1s_per_loci[nloci_after_filtering]);
            hapEntries[nloci_after_filtering].positions2.reserve(number_of_2s_per_loci[nloci_after_filtering]);
        
        }else{
            if(benchmark_flag == "XOR"){
                hapEntries[nloci_after_filtering].xors.reserve(number_of_1s_per_loci[nloci_after_filtering]);
            }
            hapEntries[nloci_after_filtering].positions.reserve(number_of_1s_per_loci[nloci_after_filtering]);
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
                    if(LOW_MEM){
                        this->hapEntries[nloci_after_filtering].xorbitset->set_bit(field);
                        this->hapEntries[nloci_after_filtering].xorbitset->num_1s++;
                    }else{
                        hapEntries[nloci_after_filtering].positions2.push_back(field);
                        //hapEntries[nloci_after_filtering].count2++;
                    }
                }
                else if (allele1 == '1' || allele2 == '1'){
                    if(LOW_MEM){
                        hapEntries[nloci_after_filtering].hapbitset->set_bit(field);
                        hapEntries[nloci_after_filtering].hapbitset->num_1s++;
                    }else{
                        hapEntries[nloci_after_filtering].positions.push_back(field);
                        //hapEntries[nloci_after_filtering].count1++;
                    }
                }
            }else{ // phased
                if(allele1 == '1'){
                    if(LOW_MEM){
                        hapEntries[nloci_after_filtering].hapbitset->set_bit(2 * field);
                        hapEntries[nloci_after_filtering].hapbitset->num_1s++;
                    }else{
                        hapEntries[nloci_after_filtering].positions.push_back(2 * field);
                    }
                }
                if(allele2 == '1'){
                    if(LOW_MEM){
                        hapEntries[nloci_after_filtering].hapbitset->set_bit(2 * field + 1);
                        hapEntries[nloci_after_filtering].hapbitset->num_1s++;
                    }else{
                        hapEntries[nloci_after_filtering].positions.push_back(2 * field + 1);
                    }
                }   
            }
        }
        nloci_after_filtering++;
    }

    // PHASE 3:  XOR
    xor_for_phased_and_unphased();
    

    // PHASE 4:  FLIP
    // for (int locus = 0; locus < nloci_after_filtering; locus++){
    //     if(hapEntries[locus].flipped){
    //         vector<int> zero_positions(nhaps - hapEntries[locus].positions.size());
    //         int j = 0;
    //         int front_one = hapEntries[locus].positions[j++];
    //         for(int i=0; i<nhaps; i++){
    //             if(i==front_one){
    //                 front_one = hapEntries[locus].positions[j++];
    //             }else{
    //                 zero_positions.push_back(i);
    //             }   
    //         }
    //         hapEntries[locus].positions = zero_positions;
    //     }
    // }
    // //PHASE 4: FLIP

    /*
    if(benchmark_flag == "XOR" || benchmark_flag == "FLIP_ONLY" ){
        for (int locus_after_filter = 0; locus_after_filter < this->nloci; locus_after_filter++){
            if(hapEntries[locus_after_filter].positions.size() > this->nhaps/2){
                hapEntries[locus_after_filter].flipped = true;

                vector<int> copy_pos;
                copy_pos.reserve(this->nhaps - hapEntries[locus_after_filter].positions.size());
                int cnt = 0;
                for(int i = 0; i< this->nhaps; i++){
                    int curr =  hapEntries[locus_after_filter].positions[cnt];
                    if(i==curr){
                        cnt++;
                    }else{
                        copy_pos.push_back(i);
                    }
                }
                
                this->hapEntries[locus_after_filter].positions = copy_pos;
                copy_pos.clear();
            }else{
                this->hapEntries[locus_after_filter].flipped = false;
            }
        }
    }
    */

    if(SKIP){
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
void HapData::readHapDataVCFMissing(string filename)
{
    float MISSING_PERCENTAGE = 0.05;
    //RELATED FLAGS: ARG_SKIP, ARG_KEEP, ARG_UNPHASED, ARG_LOW_MEM, ARG_MAF, ARG_MISSING

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
        for (int field = 0; field < current_nhaps; field++)
        {
            ss >> junk;
            allele1 = junk[0];
            separator = junk[1];
            allele2 = junk[2];
            if ( (allele1 != '0' && allele1 != '1') || (allele2 != '0' && allele2 != '1') )
            {
                // if(allele1 == '.' || allele2 == '.'){
                //     continue;
                // }
                cerr << "ERROR XX: Alleles must be coded 0/1 only.\n";
                cerr << allele1 << " " << allele2 << endl;
                //throw 0;
            }

            //if(separator != '|'){
            //    cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
            //    throw 0;
            //}

            if(unphased){
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
                if(allele1 == '.' || allele2 == '.'){
                    missing_count++;
                }
            }
        }

        int derived_allele_count = (unphased? (number_of_1s + number_of_2s*2) : number_of_1s);

        // --skip-low-freq filtering based on MAF
        if ( SKIP && (derived_allele_count*1.0/((current_nhaps-missing_count)*2) < MAF || ((current_nhaps-derived_allele_count-missing_count)*1.0/((current_nhaps-missing_count)*2)) < MAF ) || 1.0*missing_count/current_nhaps > MISSING_PERCENTAGE ) {
            skiplist.push(nloci_before_filtering-1);
            skipcount++;
        } else {
            number_of_1s_per_loci.push_back(number_of_1s);
            if(unphased){
                number_of_2s_per_loci.push_back(number_of_2s);
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

    if(SKIP){
        cerr << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
        (*flog)  << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
    }
   
    int nhaps = unphased ? (current_nhaps ) : (current_nhaps ) * 2;
    cerr << "Loading " << nhaps << " haplotypes and " << nloci_before_filtering-skipcount << " loci...\n";
    initHapData(nhaps, nloci_before_filtering-skipcount);

    string junk;
    char allele1, allele2, separator;
    bool skipLine = false; // to skip metadata lines

    skipQueue = skiplist; 
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
        
        if(unphased){
            if(benchmark_flag == "XOR"){
                hapEntries[nloci_after_filtering].xors.reserve(number_of_1s_per_loci[nloci_after_filtering]+number_of_2s_per_loci[nloci_after_filtering]);
            }

            hapEntries[nloci_after_filtering].positions.reserve(number_of_1s_per_loci[nloci_after_filtering]);
            hapEntries[nloci_after_filtering].positions2.reserve(number_of_2s_per_loci[nloci_after_filtering]);
        
        }else{
            if(benchmark_flag == "XOR"){
                hapEntries[nloci_after_filtering].xors.reserve(number_of_1s_per_loci[nloci_after_filtering]);
            }
            hapEntries[nloci_after_filtering].positions.reserve(number_of_1s_per_loci[nloci_after_filtering]);
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
                    if(LOW_MEM){
                        this->hapEntries[nloci_after_filtering].xorbitset->set_bit(field);
                        this->hapEntries[nloci_after_filtering].xorbitset->num_1s++;
                    }else{
                        hapEntries[nloci_after_filtering].positions2.push_back(field);
                        //hapEntries[nloci_after_filtering].count2++;
                    }
                }
                else if (allele1 == '1' || allele2 == '1'){
                    if(LOW_MEM){
                        hapEntries[nloci_after_filtering].hapbitset->set_bit(field);
                        hapEntries[nloci_after_filtering].hapbitset->num_1s++;
                    }else{
                        hapEntries[nloci_after_filtering].positions.push_back(field);
                        //hapEntries[nloci_after_filtering].count1++;
                    }
                }
            }else{ // phased
                if(allele1 == '1'){
                    if(LOW_MEM){
                        hapEntries[nloci_after_filtering].hapbitset->set_bit(2 * field);
                        hapEntries[nloci_after_filtering].hapbitset->num_1s++;
                    }else{
                        hapEntries[nloci_after_filtering].positions.push_back(2 * field);
                    }
                }
                if(allele2 == '1'){
                    if(LOW_MEM){
                        hapEntries[nloci_after_filtering].hapbitset->set_bit(2 * field + 1);
                        hapEntries[nloci_after_filtering].hapbitset->num_1s++;
                    }else{
                        hapEntries[nloci_after_filtering].positions.push_back(2 * field + 1);
                    }
                }   
                if(allele1 == '.' || allele2 == '.'){
                    if(LOW_MEM){
                        hapEntries[nloci_after_filtering].xorbitset->set_bit(2 * field + 1);
                        hapEntries[nloci_after_filtering].xorbitset->num_1s++;
                    }else{
                        throw 1;
                    }
                }
            }
        }
        nloci_after_filtering++;
    }

    // PHASE 3:  XOR
    if(!MISSING)
        xor_for_phased_and_unphased();
    

    // PHASE 4:  FLIP
    // for (int locus = 0; locus < nloci_after_filtering; locus++){
    //     if(hapEntries[locus].flipped){
    //         vector<int> zero_positions(nhaps - hapEntries[locus].positions.size());
    //         int j = 0;
    //         int front_one = hapEntries[locus].positions[j++];
    //         for(int i=0; i<nhaps; i++){
    //             if(i==front_one){
    //                 front_one = hapEntries[locus].positions[j++];
    //             }else{
    //                 zero_positions.push_back(i);
    //             }   
    //         }
    //         hapEntries[locus].positions = zero_positions;
    //     }
    // }
    // //PHASE 4: FLIP

    /*
    if(benchmark_flag == "XOR" || benchmark_flag == "FLIP_ONLY" ){
        for (int locus_after_filter = 0; locus_after_filter < this->nloci; locus_after_filter++){
            if(hapEntries[locus_after_filter].positions.size() > this->nhaps/2){
                hapEntries[locus_after_filter].flipped = true;

                vector<int> copy_pos;
                copy_pos.reserve(this->nhaps - hapEntries[locus_after_filter].positions.size());
                int cnt = 0;
                for(int i = 0; i< this->nhaps; i++){
                    int curr =  hapEntries[locus_after_filter].positions[cnt];
                    if(i==curr){
                        cnt++;
                    }else{
                        copy_pos.push_back(i);
                    }
                }
                
                this->hapEntries[locus_after_filter].positions = copy_pos;
                copy_pos.clear();
            }else{
                this->hapEntries[locus_after_filter].flipped = false;
            }
        }
    }
    */

    if(SKIP){
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

void HapData::readHapDataTPED(string filename)
{
    if(unphased){
        cerr<<"TPED for unphased not implemented"<<endl;
        throw 1;
    }
    
    igzstream fin;
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
    int nloci_before_filter = 0;
    int previous_nhaps = -1;
    int current_nhaps = 0;
    //Counts number of haps (cols) and number of loci (rows)
    //if any lines differ, send an error message and throw an exception
    while (getline(fin, line))
    {
        //getline(fin,line);
        //if(fin.eof()) break;
        nloci_before_filter++;
        current_nhaps = countFields(line);
        
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
    }

    fin.clear(); // clear error flags
    //fin.seekg(fileStart);
    fin.close();
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    cerr << "Loading " << current_nhaps - numMapCols << " haplotypes and " << nloci_before_filter << " loci...\n";
    
    if (unphased){
        initHapData((current_nhaps - numMapCols)/2, nloci_before_filter); //TODO
    }


    vector<vector<int> > positions(nloci_before_filter);
    string junk;
    char allele1, allele2;
    queue<int> skipQueue_local;
    for (int locus = 0; locus < nloci_before_filter; locus++)
    {
        for (int i = 0; i < numMapCols; i++)
        {
            fin >> junk;
        }
        for (int hap = 0; hap < current_nhaps - numMapCols; hap++)
        {
            if(unphased){
                fin >> allele1;
                fin >> allele2;
                if ( (allele1 != '0' && allele1 != '1') || (allele2 != '0' && allele2 != '1') ){
                    cerr << "ERROR: Alleles must be coded 0/1 only.\n";
                    cerr << allele1 << " " << allele2 << endl;
                    throw 0;
                }

                char allele = '0';
                if (allele1 == '1' && allele2 == '1'){
                    allele = '2';
                    hapEntries[locus].positions2.push_back(hap);
                    //hapEntries[locus].positions.push_back(2*hap);
                }
                else if (allele1 == '1' || allele2 == '1'){
                    allele = '1';
                    hapEntries[locus].positions.push_back(hap);
                    //hapEntries[locus].positions.push_back(2*hap+1);

                }
                //else{
                //allele = '0' implied;
                //}
                //data->data[hap][locus] = allele;
                
            }
            else{
                char allele;
                fin >> allele;
                //fin >> data->data[hap][locus];
                if (allele!= '0' && allele!= '1')
                {
                    cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
                    cerr << "Allele "<<allele << endl;
                    throw 0;
                }
                if(allele=='1'){
                    positions[locus].push_back(hap);
                }
            }
        }
        if( SKIP && (positions[locus].size()*1.0/current_nhaps < MAF || 1-(positions[locus].size()*1.0/current_nhaps) < MAF ) ) {
            skipQueue_local.push(locus);
            vector<int>().swap(positions[locus]);
        }
    }

    fin.close();

    if(!unphased){
        initHapData(current_nhaps - numMapCols, nloci_before_filter-skipQueue_local.size());
    }
    this->skipQueue =  queue<int>(skipQueue_local); 
    int loc_after_filter=0;
    for(int i=0; i<nloci_before_filter; i++){
         if(!skipQueue_local.empty()){
            if(skipQueue_local.front() == i){
                skipQueue_local.pop();    
                this->skipQueue.push(i); // make a copy to skipQueue
                continue;
            }
        }
        hapEntries[loc_after_filter++].positions = positions[i];
    }

}

