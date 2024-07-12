#include "hapdata.h"
#include "../gzstream.h"
#include "../selscan-cli.h"
#include <algorithm>
#include <sstream>

void HapData::initHapData_bitset(int nhaps, unsigned int nloci)
{
    if (nhaps < 1 || nloci < 1)
    {
        cerr << "ERROR: number of haplotypes (" << nhaps << ") and number of loci (" << nloci << ") must be positive.\n";
        throw 0;
    }

    this->hapEntries = new struct HapEntry[nloci];
    this->nhaps = nhaps;
    this->nloci = nloci;
         cout << "ERROR: number of haplotypes (" << nhaps << ") and number of loci (" << nloci << ") must be positive.\n";

    for (unsigned int j = 0; j < nloci; j++){
        hapEntries[j].hapbitset = new MyBitset(nhaps);
        hapEntries[j].xorbitset = new MyBitset(nhaps);
    }
    INIT_SUCCESS = true;
}

void HapData::releaseHapData_bitset()
{
    if (hapEntries == NULL) return;

    //we have a MyBitset for every locus
    for (unsigned int j = 0; j < nloci; j++){
        delete hapEntries[j].hapbitset ; //MyBitset destructor called
        delete hapEntries[j].xorbitset; //MyBitset destructor called
    }
    delete [] hapEntries;
    hapEntries = NULL;
    this->nhaps = -9;
    this->nloci = -9;
    return;
}




/** Sets up structure according to nhaps and nloci
 * 
*/
void HapData::initHapData(int nhaps, int nloci)
{
    if (nhaps < 1 || nloci < 1)
    {
        cerr << "ERROR: number of haplotypes (" << nhaps << ") and number of loci (" << nloci << ") must be positive.\n";
        throw 0;
    }
    cout<<"hapdata with nhaps: "<<nhaps<<"nloci "<<nloci<<endl;
    this->hapEntries = new struct HapEntry[nloci];
    this->nhaps = nhaps;
    this->nloci = nloci;
    INIT_SUCCESS = true;
}

void HapData::releaseHapData()
{
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


// START BITSET
void HapData::readHapData_bitset(string filename)
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

    //int fileStart = fin.tellg();
    string line;
    int previous_nhaps = -1;
    int current_nhaps = 0;
    
    //Counts number of haps (rows) and number of loci (cols)
    //if any lines differ, send an error message and throw an exception

    queue<int> skiplist;
    vector<int> num_1s_per_loci;
    int nloci_before_filter = 0;
    while (getline(fin, line))
    {
        //getline(fin,line);
        //if(fin.eof()) break;
        
        nloci_before_filter++;
        pair<int, int> fo = countFieldsAndOnes(line);
        current_nhaps = fo.first;
        num_1s_per_loci.push_back(fo.second);
        if( SKIP && (fo.second*1.0/current_nhaps < MAF || 1-(fo.second*1.0/current_nhaps) < MAF ) ) {
            skiplist.push(nloci_before_filter-1);
        }

        //cout << "nloci: " << current_nloci << endl;
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
    //fin.seekg(fileStart);
    fin.close();


    //PHASE 2: Open VCF File To Load into Data Structure
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    cerr << "Loading " << current_nhaps << " haplotypes and " << nloci_before_filter << " loci...\n";
    

    if (unphased){
        initHapData_bitset(current_nhaps/2, nloci_before_filter-skiplist.size());
    }
    else{
        initHapData_bitset(current_nhaps, nloci_before_filter-skiplist.size());
    }

    this->skipQueue = queue<int>();
    

    char allele1;
    //int locus_after_filter = 0;
    for (int locus = 0; locus < this->nloci; locus++)
    {
        if(!skiplist.empty()){
            if(skiplist.front() == locus){
                skiplist.pop();
                skipQueue.push(locus);
                cout<<"skipping locus "<<locus<<endl;
                getline(fin, line);
                continue;
            }
        }

        this->hapEntries[locus].hapbitset->num_1s = num_1s_per_loci[locus];

        vector<bool> current_haps(current_nhaps, false);
        for (int hap = 0; hap < current_nhaps; hap++)
        {
            if(unphased){
                fin >> allele1;
                if (allele1 != '0' && allele1 != '1'){
                    cerr << "ERROR: Alleles must be coded 0/1 only.\n";
                    cerr << allele1 << endl;
                    throw 0;
                }
                current_haps[hap] = (allele1 == '1');
            }
            else{
                fin >> allele1;
                if (allele1 != '0' && allele1 != '1')
                {
                    cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
                    throw 0;
                }
                if(allele1=='1'){
                    
                    this->hapEntries[locus].hapbitset->set_bit(hap);
                }
            }
        }

        if (unphased){
            cerr<<"Unphased lowmem not implemented"<<endl;
            throw 0;
            if (current_nhaps % 2 != 0)
            {
                cerr << "ERROR:  Number of haplotypes must be even for unphased.\n";
                throw 0;
            }

            for (int hap = 0; hap < current_nhaps/2; hap++){ 
                if (hap % 2 == 1){
                    if (current_haps[hap]  && current_haps[hap*2]){ 
                        //data->data[(hap-1)/2][locus] = '2';
                        this->hapEntries[locus].positions2.push_back(hap);
                        //this->hapEntries[locus].positions.push_back(2*hap);
                    }
                    else if ( (current_haps[hap] && !current_haps[hap*2]) || (!current_haps[hap] && current_haps[hap*2]) ){
                        //data->data[(hap-1)/2][locus] = '1';
                        this->hapEntries[locus].positions.push_back(hap);
                        //this->hapEntries[locus].positions.push_back(2*hap+1);

                    }
                }
                // else{
                //     data->data[hap/2][locus] = allele1;
                // }
            }
        }
    }
    fin.close();



    //PHASE 3: XOR
    for(int locus_after_filter = 0; locus_after_filter < this->nloci; locus_after_filter++){
        if(locus_after_filter==0){
            MyBitset* b1 =(hapEntries[locus_after_filter].hapbitset);
            for (int k = 0; k < b1->nwords; k++) {
                hapEntries[locus_after_filter].xorbitset->bits[k] = b1->bits[k] ;
            }
            hapEntries[locus_after_filter].xorbitset->num_1s = b1->num_1s;
        }else{
            MyBitset* b1 =(hapEntries[locus_after_filter].hapbitset);
            MyBitset* b2 = (hapEntries[locus_after_filter-1].hapbitset);

            int sum = 0;
            for (int k = 0; k < b1->nwords; k++) {
                hapEntries[locus_after_filter].xorbitset->bits[k] = b1->bits[k] ^ b2->bits[k];
                sum += __builtin_popcountll(hapEntries[locus_after_filter].xorbitset->bits[k]);
            }
            hapEntries[locus_after_filter].xorbitset->num_1s = sum;
            
        }
    }

    //PHASE 4: FLIP
    // for (int locus_after_filter = 0; locus_after_filter < this->nloci; locus_after_filter++){
    //     if(hapEntries[locus_after_filter].hapbitset->num_1s > nhaps/2){
    //         hapEntries[locus_after_filter].flipped = true;
    //         MyBitset* b1;
    //         b1 = hapEntries[locus_after_filter].hapbitset;

    //         //#pragma omp simd
    //         for(int k = 0; k<b1->nwords; k++){
    //             b1->bits[k] = ~(b1->bits[k]);   // negate all bits
    //         }

    //         //#pragma omp simd
    //         for(int i = b1->nbits; i<b1->nwords*b1->WORDSZ; i++){
    //             b1->clear_bit(i);       // clear the trailing bits
    //         }
    //     }else{
    //         hapEntries[locus_after_filter].flipped = false;
    //     }
    // }
    for(int locus_after_filter = 0; locus_after_filter < this->nloci; locus_after_filter++){
        cout<<locus_after_filter<<"::: ";
        hapEntries[locus_after_filter].hapbitset->print_pos();
    }
    exit(1);
}


void HapData::readHapDataVCF_bitset(string filename)
{
    igzstream fin;

    vector<int> number_of_1s_per_loci;
    vector<int> number_of_2s_per_loci;
    queue<int> skiplist;

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

    //Counts number of haps (cols) and number of loci (rows)
    //if any lines differ, send an error message and throw an exception
    int num_meta_data_lines = 0;
    while (getline(fin, line))
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
                    //allele = '2';
                    number_of_2s++;
                }
                else if (allele1 == '1' || allele2 == '1'){
                    number_of_1s++;

                }else{
                    // allele = '0' implied;
                }
                //data->data[field][locus] = allele;
            }
            else{
                if(allele1 == '1'){
                    number_of_1s++;
                }
                if(allele2 == '1'){
                    number_of_1s++;
                }
                // data->data[2 * field][locus] = allele1;
                // data->data[2 * field + 1][locus] = allele2;
            }
        }

        int derived_allele_count = (unphased? (number_of_1s + number_of_2s*2) : number_of_1s);

        if ( SKIP && (derived_allele_count*1.0/(current_nhaps*2) < MAF || 1-(derived_allele_count*1.0/(current_nhaps*2)) < MAF ) ) {
            skiplist.push(nloci_before_filtering-1);
            skipcount++;
        } else {
            number_of_1s_per_loci.push_back(number_of_1s);
            number_of_2s_per_loci.push_back(number_of_2s);
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

    //Pass 2: Load according to first pass information
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    int nhaps = unphased ? (current_nhaps ) : (current_nhaps ) * 2;

    cerr << "Loading " << nhaps << " haplotypes and " << nloci_before_filtering << " loci...\n";
    if(SKIP){
        cerr << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
        (*flog)  << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
    }

    string junk;
    char allele1, allele2, separator;
    bool skipLine = false; // to skip metadata lines

    initHapData_bitset(nhaps, nloci_before_filtering-skipcount);

    skipQueue = skiplist; 
    unsigned int nloci_after_filtering = 0;

    for (unsigned int locus = 0; locus < nloci_before_filtering; locus++)
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
            //TODO
        }else{
            hapEntries[nloci_after_filtering].hapbitset->num_1s = number_of_1s_per_loci[nloci_after_filtering];
            if(number_of_1s_per_loci[nloci_after_filtering] > nhaps/2){
                hapEntries[nloci_after_filtering].flipped = true;
            }else{
                hapEntries[nloci_after_filtering].flipped = false;
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
                throw 0;
            }

            //if(separator != '|'){
            //    cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
            //    throw 0;
            //}
            if(unphased){
                //TODO
            }
            else{
                if(allele1 == '1'){
                    (hapEntries[nloci_after_filtering].hapbitset)->set_bit(2 * field);
                }
                if(allele2 == '1'){
                    (hapEntries[nloci_after_filtering].hapbitset)->set_bit(2 * field + 1);;
                }
            }
        }

        if(nloci_after_filtering==0){
            MyBitset* b1 =(hapEntries[nloci_after_filtering].hapbitset);
            int sum = 0;
            //#pragma omp simd
            for (int k = 0; k < b1->nwords; k++) {
                hapEntries[nloci_after_filtering].xorbitset->bits[k] = b1->bits[k] ;
            }
            hapEntries[nloci_after_filtering].xorbitset->num_1s = b1->num_1s;
        }else{
            MyBitset* b1 =(hapEntries[nloci_after_filtering].hapbitset);
            MyBitset* b2 = (hapEntries[nloci_after_filtering-1].hapbitset);

            int sum = 0;
            //#pragma omp simd
            for (int k = 0; k < b1->nwords; k++) {
                hapEntries[nloci_after_filtering].xorbitset->bits[k] = b1->bits[k] ^ b2->bits[k];
                sum += __builtin_popcountll(hapEntries[nloci_after_filtering].xorbitset->bits[k]);
            }
            hapEntries[nloci_after_filtering].xorbitset->num_1s = sum;
            
        }
        nloci_after_filtering++;
    }

    //handle fliiped
    // for (int locus = 0; locus < nloci_after_filtering; locus++){
    //     if(hapEntries[locus].flipped){
    //         MyBitset* b1;
    //         b1 = hapEntries[locus].hapbitset;

    //         //#pragma omp simd
    //         for(int k = 0; k<b1->nwords; k++){
    //                 b1->bits[k] = ~(b1->bits[k]);
    //         }

    //         //#pragma omp simd
    //         for(int i = b1->nbits; i<b1->nwords*b1->WORDSZ; i++){
    //             b1->clear_bit(i);
    //         }
    //     }
    // }


    if(SKIP){
        cerr << "Removed " << skipcount << " low frequency variants.\n";
        (*flog) << "Removed " << skipcount << " low frequency variants.\n";
    }

    fin.close();
}


void HapData::do_xor(){
    for(int i = 0; i<this->nloci; i++){
        if(i==0){
            if(benchmark_flag == "XOR"){
                hapEntries[i].xors = hapEntries[i].positions;
                // vector<unsigned int>& source = hapEntries[i].positions;
                // vector<unsigned int>& destination = hapEntries[i].xors;
                // std::copy(source.begin(), source.end(), destination.begin());
                // hapEntries[i].xors1 = hapEntries[i].positions;
                // hapEntries[i].xors2 = hapEntries[i].positions2;
            }
        }else{
            if(benchmark_flag == "XOR"){
                vector<unsigned int>& curr_xor = hapEntries[i].xors;
                vector<unsigned int>& curr_positions = hapEntries[i].positions;
                vector<unsigned int>& prev_positions = hapEntries[i-1].positions;
                std::set_symmetric_difference(curr_positions.begin(), curr_positions.end(),prev_positions.begin(), prev_positions.end(),
                                std::back_inserter(curr_xor));
            }
        }
    }
}



void HapData::readHapDataVCF(string filename)
{
    igzstream fin;
    vector<int> number_of_1s_per_loci;
    vector<int> number_of_2s_per_loci;
    queue<int> skiplist;

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
                    //allele = '2';
                    number_of_2s++;
                }
                else if (allele1 == '1' || allele2 == '1'){
                    number_of_1s++;

                }else{
                    // allele = '0' implied;
                }
                //data->data[field][locus] = allele;
            }
            else{
                if(allele1 == '1'){
                    number_of_1s++;
                }
                if(allele2 == '1'){
                    number_of_1s++;
                }
                // data->data[2 * field][locus] = allele1;
                // data->data[2 * field + 1][locus] = allele2;
            }
        }

        int derived_allele_count = (unphased? (number_of_1s + number_of_2s*2) : number_of_1s);

        if ( SKIP && (derived_allele_count*1.0/(current_nhaps*2) < MAF || 1-(derived_allele_count*1.0/(current_nhaps*2)) < MAF ) ) {
            skiplist.push(nloci_before_filtering-1);
            skipcount++;
        } else {
            number_of_1s_per_loci.push_back(number_of_1s);
            number_of_2s_per_loci.push_back(number_of_2s);
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

    //Pass 2: Load according to first pass information
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    int nhaps = unphased ? (current_nhaps ) : (current_nhaps ) * 2;

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


        // PHASE 3:  XOR
        if(unphased){
            if(nloci_after_filtering==0){
                //hapEntries[nloci_after_filtering].xors;
                //hapEntries[nloci_after_filtering].xors1 = hapEntries[nloci_after_filtering].positions;
                //hapEntries[nloci_after_filtering].xors2 = hapEntries[nloci_after_filtering].positions2;
            }else{
                getThreeUnphasedGroups(hapEntries[nloci_after_filtering].positions, hapEntries[nloci_after_filtering-1].positions,hapEntries[nloci_after_filtering].positions2, hapEntries[nloci_after_filtering-1].positions2, hapEntries[nloci_after_filtering].g);
            }

        }else{
            if(nloci_after_filtering==0){
                if(benchmark_flag == "XOR"){
                    vector<unsigned int>& source = hapEntries[nloci_after_filtering].positions;
                    vector<unsigned int>& destination = hapEntries[nloci_after_filtering].xors;
                    std::copy(source.begin(), source.end(), destination.begin());
                    hapEntries[nloci_after_filtering].xors1 = hapEntries[nloci_after_filtering].positions;
                    hapEntries[nloci_after_filtering].xors2 = hapEntries[nloci_after_filtering].positions2;
                }
            }else{
                if(benchmark_flag == "XOR"){
                    vector<unsigned int>& curr_xor = hapEntries[nloci_after_filtering].xors;
                    vector<unsigned int>& curr_positions = hapEntries[nloci_after_filtering].positions;
                    vector<unsigned int>& prev_positions = hapEntries[nloci_after_filtering-1].positions;
                    std::set_symmetric_difference(curr_positions.begin(), curr_positions.end(),prev_positions.begin(), prev_positions.end(),
                                    std::back_inserter(curr_xor));
                
                    //unphased
                    if(unphased){
                        for(int i=0; i<curr_loc_str.length(); i++){
                            int sub = curr_loc_str[i]-prev_loc_str[i];
                            if(sub<0){
                                sub = 3+sub;
                            }
                            if(sub==1){
                                hapEntries[nloci_after_filtering].xors1.push_back(i);
                            }else if(sub==2){
                                hapEntries[nloci_after_filtering].xors2.push_back(i);
                            }
                        }
                    }   
                }
            }

        }
        
        prev_loc_str = curr_loc_str;
        nloci_after_filtering++;
    }


    // PHASE 4:  FLIP
    // for (int locus = 0; locus < nloci_after_filtering; locus++){
    //     if(hapEntries[locus].flipped){
    //         vector<unsigned int> zero_positions(nhaps - hapEntries[locus].positions.size());
    //         int j = 0;
    //         unsigned int front_one = hapEntries[locus].positions[j++];
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

                vector<unsigned int> copy_pos;
                copy_pos.reserve(this->nhaps - hapEntries[locus_after_filter].positions.size());
                int cnt = 0;
                for(int i = 0; i< this->nhaps; i++){
                    unsigned int curr =  hapEntries[locus_after_filter].positions[cnt];
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
    unsigned int nloci_before_filter = 0;
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


    vector<vector<unsigned int> > positions(nloci_before_filter);
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
            vector<unsigned int>().swap(positions[locus]);
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

//reads in haplotype data and also does basic checks on integrity of format
//returns a populated HaplotypeData structure if successful
//throws an exception otherwise
//impute hap IMPUTE HAP is transposed format (thap) where row represents loci,  column replesent individual
//so wc -l of impute hap is same as map.

void HapData::readHapData(string filename)
{
    //PHASE 1: Read Hap File to get "nloci", "nhaps" and "skiplist"
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
    vector<int> num_1s_per_loci;
    int nloci_before_filter = 0;

    while (getline(fin, line)) //Counts number of haps (rows) and number of loci (cols) 
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
        else if (previous_nhaps != current_nhaps) //if any lines differ, send an error message and throw an exception
        {
            cerr << "ERROR: line " << nloci_before_filter << " of " << filename << " has " << current_nhaps
                 << ", but the previous line has " << previous_nhaps << ".\n";
            throw 0;
        }
        previous_nhaps = current_nhaps;
    }
    fin.clear();
    fin.close();


    //PHASE 2: Open Hap File To Load into Data Structure
    fin.open(filename.c_str());
    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }
    cerr << "Loading " << current_nhaps << " haplotypes and " << nloci_before_filter << " loci...\n";
    
    if (unphased){
        initHapData(current_nhaps/2, nloci_before_filter-skiplist.size());
    }
    else{
        initHapData(current_nhaps, nloci_before_filter-skiplist.size());
    }
    this->skipQueue = queue<int>(); // make a copy
    

    char allele1;

    int locus_after_filter = 0;
    getline(fin, line);
    for (int locus_before_filter = 0; locus_before_filter < nloci_before_filter; locus_before_filter++)
    {
        if(!skiplist.empty()){
            if(skiplist.front() == locus_before_filter){
                skiplist.pop();
                 this->skipQueue.push(locus_before_filter);
                getline(fin, line);
                continue;
            }
        }
        
        stringstream ss(line);
        this->hapEntries[locus_after_filter].positions.reserve(num_1s_per_loci[locus_before_filter]);
        for (int hap = 0; hap < current_nhaps; hap++)
        {
            if(unphased){
                cerr << "ERROR: UNPHASED HAP NOT IMPLEMENTED.\n";
                throw 0;

                ss >> allele1;
                if (allele1 != '0' && allele1 != '1'){
                    cerr << "ERROR: Alleles must be coded 0/1 only.\n";
                    cerr << allele1 << endl;
                    throw 0;
                }
            }
            else{
                ss >> allele1;
                if (allele1 != '0' && allele1 != '1')
                {
                    cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
                    throw 0;
                }
                if(allele1=='1'){
                    this->hapEntries[locus_after_filter].positions.push_back(hap);
                }
            }
        }
        locus_after_filter++;
        getline(fin, line);
    }
    fin.clear();
    fin.close();


    //PHASE 3: XOR

    hapEntries[0].xors = hapEntries[0].positions;
    for(int locus_after_filter = 1; locus_after_filter < this->nloci; locus_after_filter++){
        vector<unsigned int>& curr_xor = hapEntries[locus_after_filter].xors;
        vector<unsigned int>& curr_positions = hapEntries[locus_after_filter].positions;
        vector<unsigned int>& prev_positions = hapEntries[locus_after_filter-1].positions;
        std::set_symmetric_difference(curr_positions.begin(), curr_positions.end(),prev_positions.begin(), prev_positions.end(), std::back_inserter(curr_xor));  
    }

    // //PHASE 4: FLIP
    // for (int locus_after_filter = 0; locus_after_filter < this->nloci; locus_after_filter++){
    //     if(hapEntries[locus_after_filter].positions.size() > this->nhaps/2){
    //         hapEntries[locus_after_filter].flipped = true;

    //         vector<unsigned int> copy_pos;
    //         copy_pos.reserve(this->nhaps - hapEntries[locus_after_filter].positions.size());
    //         int cnt = 0;
    //         for(int i = 0; i< this->nhaps; i++){
    //             unsigned int curr =  hapEntries[locus_after_filter].positions[cnt];
    //             if(i==curr){
    //                 cnt++;
    //             }else{
    //                 copy_pos.push_back(i);
    //             }
    //         }
            
    //         this->hapEntries[locus_after_filter].positions = copy_pos;
    //         copy_pos.clear();
    //         // vector<unsigned int> zero_positions(this->nhaps - this->hapEntries[locus_after_filter].positions.size());
    //         // int j = 0;
    //         // unsigned int front_one = this->hapEntries[locus_after_filter].positions[j++];
    //         // for(int i=0; i<nhaps; i++){
    //         //     if(i==front_one){
    //         //         front_one = this->hapEntries[locus_after_filter].positions[j++];
    //         //     }else{
    //         //         zero_positions.push_back(i);
    //         //     }   
    //         // }
    //         // this->hapEntries[locus_after_filter].positions = zero_positions;
    //     }else{
    //         this->hapEntries[locus_after_filter].flipped = false;
    //     }
    // }

}




