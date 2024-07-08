#ifndef __DATA_READER_H__
#define __DATA_READER_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cctype>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <cstdio>
#include "../selscan-cli.h"
#include "../hapmap/hapdata.h"

using namespace std;

//class DataReaderUnphased
// if u get the vectors at first

class DataReader{
    public:
        const string filename;
        HapData& hapData;
        ofstream* flog;
        int num_threads;

        DataReader(const std::string& filename, HapData& hapData)
            : filename(filename), hapData(hapData), flog(hapData.flog), num_threads(hapData.num_threads)
        {}

        /// Assigns value of HapData.nhaps and HapData.nloci  
        virtual void init_based_on_lines() {}
 
        /// @brief Assigns value of #hapData.positions  #hapData.positions2  #skiplist  #hapData.nhaps
        virtual void populate_positions_skipqueue() {}
        virtual void do_xor() {}

        virtual ~DataReader() {}
};



/*
    // //now you can call init hap data?
    // void initialHapDataReading(){
    //     //read the input
    //     // get info about
    //     // nloci, nhaps, queue, num1s, num2s
    //     // line id ()
    // }
    void init_based_on_lines(){ //can do parallely in vcf
       
        igzstream fin;

        //return stuff
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

                if(unphased){
                    char allele = '0';
                    if (allele1 == '1' && allele2 == '1'){ //allele = '2';
                        number_of_2s++;
                    }
                    else if (allele1 == '1' || allele2 == '1'){
                        number_of_1s++;
                    }
                }
                else{
                    if(allele1 == '1'){
                        number_of_1s++;
                    }
                    if(allele2 == '1'){
                        number_of_1s++;
                    }
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
        fin.clear(); 
        fin.close();
        hapData.nloci = nloci_before_filtering;
        hapData.nhaps = unphased ? (current_nhaps ) : (current_nhaps ) * 2;
    }
    void n1s_n2s(){ //can do parallely in vcf

    }
    void populate_positions_skipqueue(){ //can do parallely in vcf
        //Pass 2: Load according to first pass information
        igzstream fin;
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
    void init(){ //serial

    }
    void loadGenotypes(){ //can do parallely in vcf

    }
    void doXor(){ // can do parallely

    }
    //if skip low freq populate_positions_skipqueue else n1s_n2s
    //then init


    void readHapDataSerial(string filename, HapData & hapData )
    {
        

        


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

        if(SKIP){
            cerr << "Removed " << skipcount << " low frequency variants.\n";
            (*flog) << "Removed " << skipcount << " low frequency variants.\n";
        }

        fin.close();

    }
*/

// class ParallelDataReader : public DataReader {  // not gzipped
//     void initialHapDataReading();
// };
// class SerialDataReader : public DataReader{  // gzipped
//     void initialHapDataReading();
// };



#endif