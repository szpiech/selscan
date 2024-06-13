/* selscan -- a program to calculate EHH-based scans for positive selection in genomes
   Copyright (C) 2014  Zachary A Szpiech

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
*/
#include "selscan-data.h"
#include "selscan-cli.h"
#include <sstream>
#include <queue> 
#include<algorithm>
using namespace std;

//Returns 1 on error
bool HapMap::loadHapMapData(param_main &p, int argc, char *argv[], ofstream* flog, ofstream* fout){
    this->flog = flog;
    hapData.flog = flog;

    this->fout = fout;

    string hapFilename = p.hapFilename;
    string hapFilename2 =  p.hapFilename2;
    string mapFilename = p.mapFilename;
    string tpedFilename = p.tpedFilename;
    string tpedFilename2 = p.tpedFilename2;
    string vcfFilename = p.vcfFilename;
    string vcfFilename2 = p.vcfFilename2;
    
    bool TPED = false;
    if (tpedFilename.compare(DEFAULT_FILENAME_POP1_TPED) != 0) TPED = true;

    bool VCF = false;
    if (vcfFilename.compare(DEFAULT_FILENAME_POP1_VCF) != 0) VCF = true;

    bool UNPHASED = p.UNPHASED;
    
    bool USE_PMAP = p.USE_PMAP;
    bool ALT = p.ALT;
    bool WAGH = p.WAGH;
    bool CALC_IHS = p.CALC_IHS;
    bool CALC_XPNSL = p.CALC_XPNSL;
    bool CALC_NSL = p.CALC_NSL;
    bool WRITE_DETAILED_IHS = p.WRITE_DETAILED_IHS;
    bool CALC_XP = p.CALC_XP;
    bool CALC_SOFT = p.CALC_SOFT;
    bool SINGLE_EHH = p.SINGLE_EHH;

    
    hapData.initParams(UNPHASED, p.SKIP, p.MAF);
    hapData2.initParams(UNPHASED, p.SKIP, p.MAF);

    try
    {
        if (TPED)
        {
            hapData.readHapDataTPED(tpedFilename);
            if (CALC_XP || CALC_XPNSL)
            {
                hapData2.readHapDataTPED(tpedFilename2);
                if (hapData.nloci != hapData2.nloci)
                {
                    std::cerr << "ERROR: Haplotypes from " << tpedFilename << " and " << tpedFilename2 << " do not have the same number of loci.\n";
                    return 1;
                }
            }
            mapData.readMapDataTPED(tpedFilename, hapData.nloci, hapData.nhaps, USE_PMAP);
        }
        else if (VCF) {
            
            hapData.readHapDataVCF(vcfFilename);
            if (CALC_XP || CALC_XPNSL)
            {
                hapData2.readHapDataVCF(vcfFilename2);
                if (hapData.nloci != hapData2.nloci)
                {
                    std::cerr << "ERROR: Haplotypes from " << vcfFilename << " and " << vcfFilename2 << " do not have the same number of loci.\n";
                    return 1;
                }
            }
            if(!CALC_NSL && !CALC_XPNSL && !USE_PMAP) {
                mapData.readMapData(mapFilename, hapData.nloci, USE_PMAP);
            }
            else{//Load physical positions
                mapData.readMapDataVCF(vcfFilename, hapData.nloci, hapData.skipQueue);
            }
        }
        else
        {
            /*
            hapData.readHapData(hapFilename);
            if (CALC_XP || CALC_XPNSL)
            {
                hapData2.readHapData(hapFilename2);
                if (hapData.nloci != hapData2.nloci)
                {
                    std::cerr << "ERROR: Haplotypes from " << hapFilename << " and " << hapFilename2 << " do not have the same number of loci.\n";
                    return 1;
                }
            }
            mapData.readMapData(mapFilename, hapData.nloci, USE_PMAP);
            */
        }
    }
    catch (...)
    {
        return 1;
    }

    //mapData.print();


    // Check if map is in order
    for (int i = 1; i < mapData.nloci; i++) {
        if ( mapData.mapEntries[i].physicalPos <  mapData.mapEntries[i-1].physicalPos ) {
            std::cerr << "ERROR: Variant physical position must be monotonically increasing.\n";
            std::cerr << "\t" << mapData.mapEntries[i].locusName << " " << mapData.mapEntries[i].physicalPos << " appears after";
            std::cerr << "\t" <<  mapData.mapEntries[i-1].locusName << " " << mapData.mapEntries[i-1].physicalPos << "\n";
            return 1;
        }
        if ( !CALC_NSL && mapData.mapEntries[i].geneticPos  < mapData.mapEntries[i-1].geneticPos  ) {
            std::cerr << "ERROR: Variant genetic position must be monotonically increasing.\n";
            std::cerr << "\t" << mapData.mapEntries[i].locusName << " " << mapData.mapEntries[i].geneticPos << " appears after";
            std::cerr << "\t" << mapData.mapEntries[i-1].locusName << " " << mapData.mapEntries[i-1].geneticPos << "\n";
            return 1;
        }
    }

    return 0;
}

/**
 * reads in map data and also does basic checks on integrity of format
 * @returns a populated MapData structure if successful
 * @throws an exception if not successful
 */
void MapData::readMapData(string filename, int expected_loci, bool USE_PMAP)
{
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
    int nloci = 0;
    int num_cols = 4;
    int current_cols = 0;
    while (getline(fin, line))
    {
        nloci++;
        current_cols = countFields(line);
        if (current_cols != num_cols)
        {
            cerr << "ERROR: line " << nloci << " of " << filename << " has " << current_cols
                 << ", but expected " << num_cols << ".\n";
            throw 0;
        }
    }

    if (nloci != expected_loci)
    {
        cerr << "ERROR: Expected " << expected_loci << " loci in map file but found " << nloci << ".\n";
        throw 0;
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

    cerr << "Loading map data for " << nloci << " loci\n";

    initMapData(nloci);

    double Mb = 1000000.0;

    string chr;
    for (unsigned int locus = 0; locus < this->nloci; locus++)
    {
        fin >> mapEntries[locus].chr;
        fin >> mapEntries[locus].locusName;
        fin >> mapEntries[locus].geneticPos;
        fin >> mapEntries[locus].physicalPos;

        locus_query_map[mapEntries[locus].locusName] = locus;


        if (USE_PMAP) mapEntries[locus].geneticPos = double(mapEntries[locus].physicalPos)/Mb;
    }

    fin.close();
}

void MapData::readMapDataTPED(string filename, int expected_loci, int expected_haps, bool USE_PMAP)
{
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
    int nloci = 0;
    int num_cols = 4;
    int current_cols = 0;
    while (getline(fin, line))
    {
        nloci++;
        current_cols = countFields(line);
        if (current_cols != num_cols + expected_haps)
        {
            cerr << "ERROR: line " << nloci << " of " << filename << " has " << current_cols
                 << ", but expected " << num_cols + expected_haps << ".\n";
            throw 0;
        }
    }

    if (nloci != expected_loci)
    {
        cerr << "ERROR: Expected " << expected_loci << " loci in map file but found " << nloci << ".\n";
        throw 0;
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

    cerr << "Loading map data for " << nloci << " loci\n";

    initMapData(nloci);
    
    double Mb = 1000000.0;
    
    string chr;
    for (int locus = 0; locus < this->nloci; locus++)
    {
        fin >> mapEntries[locus].chr;
        fin >> mapEntries[locus].locusName;
        fin >> mapEntries[locus].geneticPos;
        fin >> mapEntries[locus].physicalPos;

        locus_query_map[mapEntries[locus].locusName] = locus;

        if (USE_PMAP) mapEntries[locus].geneticPos = double(mapEntries[locus].physicalPos)/Mb;
        getline(fin, line);
    }

    fin.close();
}

void MapData::readMapDataVCF(string filename, int expected_loci, queue <int>& skipQueue) {
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    string line;
    int nloci = 0;
    int numCommentedLines = 0;
    while (getline(fin, line))
    {
        if (line[0] == '#') {
            numCommentedLines++;
        }
        else {
            nloci++;
        }
    }

    if (nloci-skipQueue.size() != expected_loci)
    {
        cerr << "ERROR: Expected " << expected_loci << " loci in file but found " << nloci-skipQueue.size() << ".\n";
        throw 0;
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

    cerr << "Loading map data for " << nloci-skipQueue.size() << " loci\n";

    for (int i = 0; i < numCommentedLines; i++) {
        getline(fin, line);
    }

    this->initMapData(nloci-skipQueue.size()); 
    
    double Mb = 1000000.0;
    
    string chr;

    int locus = 0;

    for (int i = 0; i < nloci; i++)
    {
        if(!skipQueue.empty()){
            if(skipQueue.front()==locus){
                skipQueue.pop();
                string junk;
                fin >> junk;
                fin >>  junk;
                fin >>  junk;
                getline(fin, line);
                continue;
            }
        }
        
        
        fin >> mapEntries[locus].chr;
        fin >> mapEntries[locus].physicalPos;
        fin >> mapEntries[locus].locusName;

        //locus_query_map[mapEntries[locus].locusName] = locus;
        locus_query_map[to_string(mapEntries[locus].physicalPos)] = locus;
        mapEntries[locus].locId = i;

        getline(fin, line);
        mapEntries[locus].geneticPos = double(mapEntries[locus].physicalPos)/Mb;

        locus++;
    }

    fin.close();
}

//allocates the arrays and populates them with MISSING or "--" depending on type
void MapData::initMapData(int nloci)
{
    if (nloci < 1)
    {
        cerr << "ERROR: number of loci (" << nloci << ") must be positive.\n";
        throw 0;
    }

    //MapData *data = new MapData;
    //

    mapEntries = new struct MapEntry[nloci];
    this->nloci = nloci;


    // data->nloci = nloci;
    // data->locusName = new string[nloci];
    // data->physicalPos = new int[nloci];
    // data->geneticPos = new double[nloci];

    // for (int locus = 0; locus < nloci; locus++)
    // {
    //     data->locusName[locus] = "--";
    //     data->physicalPos[locus] = MISSING;
    //     data->geneticPos[locus] = MISSING;
    // }
}

void MapData::releaseMapData()
{
    if (mapEntries == NULL) return;
    this->nloci = -9;
    // delete [] data->locusName;
    // delete [] data->physicalPos;
    // delete [] data->geneticPos;
    // delete data;
    // data = NULL;
    
    delete [] mapEntries;
    mapEntries = NULL;
    return;
}

//reads in haplotype data and also does basic checks on integrity of format
//returns a populated HaplotypeData structure if successful
//throws an exception otherwise
//impute hap IMPUTE HAP is transposed format (thap) where row represents loci,  column replesent individual
//so wc -l of impute hap is same as map.

void HapData::readHapData(string filename)
{
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
    unsigned int nloci = 0;
    int previous_nhaps = -1;
    int current_nhaps = 0;
    //Counts number of haps (rows) and number of loci (cols)
    //if any lines differ, send an error message and throw an exception
    while (getline(fin, line))
    {
        //getline(fin,line);
        //if(fin.eof()) break;
        nloci++;

        current_nhaps = countFields(line);
        //cout << "nloci: " << current_nloci << endl;
        if (previous_nhaps < 0)
        {
            previous_nhaps = current_nhaps;
            continue;
        }
        else if (previous_nhaps != current_nhaps)
        {
            cerr << "ERROR: line " << nloci << " of " << filename << " has " << current_nhaps
                 << ", but the previous line has " << previous_nhaps << ".\n";
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

    cerr << "Loading " << current_nhaps << " haplotypes and " << nloci << " loci...\n";
    

    if (unphased){
        initHapData(current_nhaps/2, nloci);
    }
    else{
        initHapData(current_nhaps, nloci);
    }


    char allele1;
    for (int locus = 0; locus < this->nloci; locus++)
    {
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
                    this->hapEntries[locus].positions.push_back(hap);
                }
            }


        }
        if(unphased){
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
}


void HapData::readHapDataTPED(string filename)
{
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
    unsigned int nloci = 0;
    int previous_nhaps = -1;
    int current_nhaps = 0;
    //Counts number of haps (cols) and number of loci (rows)
    //if any lines differ, send an error message and throw an exception
    while (getline(fin, line))
    {
        //getline(fin,line);
        //if(fin.eof()) break;
        nloci++;
        current_nhaps = countFields(line);
        //cout << "nloci: " << current_nhaps << endl;
        if (previous_nhaps < 0)
        {
            previous_nhaps = current_nhaps;
            continue;
        }
        else if (previous_nhaps != current_nhaps)
        {
            cerr << "ERROR: line " << nloci << " of " << filename << " has " << current_nhaps
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

    cerr << "Loading " << current_nhaps - numMapCols << " haplotypes and " << nloci << " loci...\n";
    

    if (unphased){
        initHapData((current_nhaps - numMapCols)/2, nloci);
    }
    else{
        initHapData(current_nhaps - numMapCols, nloci);
    }

    string junk;
    char allele1, allele2;

    for (int locus = 0; locus < this->nloci; locus++)
    {
        for (int i = 0; i < numMapCols; i++)
        {
            fin >> junk;
        }
        for (int hap = 0; hap < this->nhaps; hap++)
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
                    throw 0;
                }
                if(allele=='1'){
                    hapEntries[locus].positions.push_back(hap);
                }
            }
        }
    }

    fin.close();
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

    initHapData(nhaps, nloci_before_filtering-skipcount);

    skipQueue = skiplist; 
    unsigned int nloci_after_filtering = 0;

    string prev_loc_str = "";
    string curr_loc_str = "";
    for (unsigned int locus = 0; locus < nloci_before_filtering+num_meta_data_lines; locus++)
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
            //locus--;
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

            if(benchmark_flag != "BITSET"){
                hapEntries[nloci_after_filtering].positions.reserve(number_of_1s_per_loci[nloci_after_filtering]);
                hapEntries[nloci_after_filtering].positions2.reserve(number_of_2s_per_loci[nloci_after_filtering]);
            }
        }else{
            if(benchmark_flag == "XOR"){
                hapEntries[nloci_after_filtering].xors.reserve(number_of_1s_per_loci[nloci_after_filtering]);
            }

            if(benchmark_flag != "BITSET"){
                hapEntries[nloci_after_filtering].positions.reserve(number_of_1s_per_loci[nloci_after_filtering]);
            }
            if(benchmark_flag == "BITSET"){
                hapEntries[nloci_after_filtering].hapbitset->num_1s = number_of_1s_per_loci[nloci_after_filtering];
                if(number_of_1s_per_loci[nloci_after_filtering] > nhaps/2){
                    hapEntries[nloci_after_filtering].flipped = true;
                }else{
                    hapEntries[nloci_after_filtering].flipped = false;
                }
            }
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
                //data->data[field][locus] = allele;
            }
            else{
                if(benchmark_flag == "BITSET"){
                    if(allele1 == '1'){
                        (hapEntries[nloci_after_filtering].hapbitset)->set_bit(2 * field);
                    }
                    if(allele2 == '1'){
                        (hapEntries[nloci_after_filtering].hapbitset)->set_bit(2 * field + 1);;
                    }
                }

                if(benchmark_flag != "BITSET"){
                    if (hapEntries[nloci_after_filtering].flipped){
                        if(allele1 == '0'){
                            hapEntries[nloci_after_filtering].positions.push_back(2 * field);
                        }
                        if(allele2 == '0'){
                            hapEntries[nloci_after_filtering].positions.push_back(2 * field + 1);
                        }
                    }else{
                        if(allele1 == '1'){
                        hapEntries[nloci_after_filtering].positions.push_back(2 * field);
                        }
                        if(allele2 == '1'){
                            hapEntries[nloci_after_filtering].positions.push_back(2 * field + 1);
                        }
                    }
                }
                
                
                // data->data[2 * field][locus] = allele1;
                // data->data[2 * field + 1][locus] = allele2;
            }
        }


        if(nloci_after_filtering==0){
            if(benchmark_flag == "XOR"){
                vector<unsigned int>& source = hapEntries[nloci_after_filtering].positions;
                vector<unsigned int>& destination = hapEntries[nloci_after_filtering].xors;
                std::copy(source.begin(), source.end(), destination.begin());

                //std::copy(hapEntries[nloci_after_filtering].positions.begin(), hapEntries[nloci_after_filtering].positions.end(), hapEntries[nloci_after_filtering].xors1.begin());
                //std::copy(hapEntries[nloci_after_filtering].positions2.begin(), hapEntries[nloci_after_filtering].positions2.end(), hapEntries[nloci_after_filtering].xors2.begin());

                hapEntries[nloci_after_filtering].xors1 = hapEntries[nloci_after_filtering].positions;
                hapEntries[nloci_after_filtering].xors2 = hapEntries[nloci_after_filtering].positions2;


            }

            if(benchmark_flag == "BITSET"){
                MyBitset* b1 =(hapEntries[nloci_after_filtering].hapbitset);
                int sum = 0;
                for (int k = 0; k < b1->nwords; k++) {
                    hapEntries[nloci_after_filtering].xorbitset->bits[k] = b1->bits[k] ;
                }
                hapEntries[nloci_after_filtering].xorbitset->num_1s = b1->num_1s;
            }
           
           
        }else{
            if(benchmark_flag == "BITSET"){
                MyBitset* b1 =(hapEntries[nloci_after_filtering].hapbitset);
                MyBitset* b2 = (hapEntries[nloci_after_filtering-1].hapbitset);

                int sum = 0;
                for (int k = 0; k < b1->nwords; k++) {
                    hapEntries[nloci_after_filtering].xorbitset->bits[k] = b1->bits[k] ^ b2->bits[k];
                    sum += __builtin_popcountll(hapEntries[nloci_after_filtering].xorbitset->bits[k]);
                }
                hapEntries[nloci_after_filtering].xorbitset->num_1s = sum;
            }
            if(benchmark_flag == "XOR"){
                vector<unsigned int>& curr_xor = hapEntries[nloci_after_filtering].xors;
                vector<unsigned int>& curr_positions = hapEntries[nloci_after_filtering].positions;
                vector<unsigned int>& prev_positions = hapEntries[nloci_after_filtering-1].positions;
                std::set_symmetric_difference(curr_positions.begin(), curr_positions.end(),prev_positions.begin(), prev_positions.end(),
                                std::back_inserter(curr_xor));
            

                            
                //unphased
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
        prev_loc_str = curr_loc_str;
        nloci_after_filtering++;
    }

    if(benchmark_flag!="BITSET"){
        for (int locus = 0; locus < nloci_before_filtering; locus++){
            if(hapEntries[locus].flipped){
                vector<unsigned int> zero_positions(nhaps - hapEntries[locus].positions.size());

                int j = 0;
                unsigned int front_one = hapEntries[locus].positions[j++];
                for(int i=0; i<nhaps; i++){
                    if(i==front_one){
                        front_one = hapEntries[locus].positions[j++];
                    }else{
                        zero_positions.push_back(i);
                    }   
                }
                hapEntries[locus].positions = zero_positions;
            }
        }
    }

    if(benchmark_flag=="BITSET"){
        for (int locus = 0; locus < nloci_before_filtering; locus++){
            if(hapEntries[locus].flipped){
                MyBitset* b1 = hapEntries[locus].hapbitset;
                for(int k = 0; k<b1->nwords; k++){
                     b1->bits[k] = ~(b1->bits[k]);
                }
                for(int i = b1->nbits; i<b1->nwords*b1->WORDSZ; i++){
                    b1->clear_bit(i);
                }
               
            }
        }
    }
    
    
    //nloci = nloci_before_filtering - skipcount*int(SKIP);
   // nloci = nloci_after_filtering;

   /*
     for (int locus = 0; locus < 1; locus++){
        hapEntries[locus].hapbitset->print_pos();
        // for (int j = 0; j < 10; j++){
        //     cout<<hapEntries[locus].positions[j]<<" ";
        // }
        // cout<<endl;
     }
    exit(1);
    */


    /*
    //iterate over all_positions vector
    int newlocus = 0;
    for (int locus = 0; locus < nloci_before_filtering; locus++)
    {
        hapEntries[locus].positions.reserve(all_positions[locus].size());
        hapEntries[locus].xors.reserve(all_positions[locus].size());

        double AF1 = all_positions[locus].size()*1.0/nhaps;
        if(unphased){
            AF1 = (all_positions[locus].size() + all_positions2[locus].size()*2.0)/nhaps;
        }

        //if(!SKIP || min(AF1, 1-AF1) >= MAF){
            for (unsigned int pos: all_positions[locus])
            {
                hapEntries[newlocus].positions.push_back(pos);
                hapEntries[newlocus].xors.push_back(pos);
                ++newlocus;
            }
            if(unphased){
                for (unsigned int pos: all_positions2[locus])
                {
                    hapEntries[newlocus].positions2.push_back(pos);
                }
            }
        //}


        if(nloci!=newlocus){
            cerr << "DEBUG ERROR: nloci != newlocus\n"<< nloci<<" "<<newlocus<<endl;
            throw 0;
        }
        
    }
    */



   //BEGIN TESTING CODES

   //END TESTING CODES

    if(SKIP){
        cerr << "Removed " << skipcount << " low frequency variants.\n";
        (*flog) << "Removed " << skipcount << " low frequency variants.\n";
    }

    fin.close();

}

/** Sets up structure according to nhaps and nloci
 * 
*/
void HapData::initHapData(int nhaps, unsigned int nloci)
{
    if (nhaps < 1 || nloci < 1)
    {
        cerr << "ERROR: number of haplotypes (" << nhaps << ") and number of loci (" << nloci << ") must be positive.\n";
        throw 0;
    }

    this->hapEntries = new struct HapEntry[nloci];
    this->nhaps = nhaps;
    this->nloci = nloci;

    if(benchmark_flag == "BITSET"){
        for (unsigned int j = 0; j < nloci; j++){
            hapEntries[j].hapbitset = new MyBitset(nhaps);
            hapEntries[j].xorbitset = new MyBitset(nhaps);

        }
    }

    //this->data = new char *[nhaps];
    // for (unsigned int i = 0; i < nhaps; i++)
    // {
    //     data->data[i] = new char[nloci];
    //     for (unsigned int j = 0; j < nloci; j++)
    //     {
    //         data->data[i][j] = MISSING_CHAR;
    //     }
    // }

    
    // for (unsigned int j = 0; j < nloci; j++)
    // {
    //     //hapEntries[j].xors.reserve(nhaps);
    //     //hapEntries[j].positions.reserve(nhaps);
    // }
    
}

void HapData::releaseHapData()
{
    if (hapEntries == NULL) return;
    for (int i = 0; i < this->nhaps; i++)
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


int countFields(const string &str)
{
    string::const_iterator it;
    int result;
    int numFields = 0;
    int seenChar = 0;
    for (it = str.begin() ; it < str.end(); it++)
    {
        result = isspace(*it);
        if (result == 0 && seenChar == 0)
        {
            numFields++;
            seenChar = 1;
        }
        else if (result != 0)
        {
            seenChar = 0;
        }
    }
    return numFields;
}

