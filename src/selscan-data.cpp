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

using namespace std;

//Returns 1 on error
bool HapMap::loadHapMapData(param_t &params, int argc, char *argv[], ofstream* flog){
    this->flog = flog;
    hapData.flog = flog;

    //this->fout = fout;

    string hapFilename = params.getStringFlag(ARG_FILENAME_POP1);
    string hapFilename2 = params.getStringFlag(ARG_FILENAME_POP2);
    string mapFilename = params.getStringFlag(ARG_FILENAME_MAP);
    string tpedFilename = params.getStringFlag(ARG_FILENAME_POP1_TPED);
    string tpedFilename2 = params.getStringFlag(ARG_FILENAME_POP2_TPED);
    string vcfFilename = params.getStringFlag(ARG_FILENAME_POP1_VCF);
    string vcfFilename2 = params.getStringFlag(ARG_FILENAME_POP2_VCF);
    
    bool TPED = false;
    if (tpedFilename.compare(DEFAULT_FILENAME_POP1_TPED) != 0) TPED = true;

    bool VCF = false;
    if (vcfFilename.compare(DEFAULT_FILENAME_POP1_VCF) != 0) VCF = true;

    bool UNPHASED = params.getBoolFlag(ARG_UNPHASED);
    hapData.unphased = UNPHASED;
    hapData2.unphased = UNPHASED;
     

    bool USE_PMAP = params.getBoolFlag(ARG_PMAP);
    bool ALT = params.getBoolFlag(ARG_ALT);
    bool WAGH = params.getBoolFlag(ARG_WAGH);
    bool CALC_IHS = params.getBoolFlag(ARG_IHS);
    bool CALC_XPNSL = params.getBoolFlag(ARG_XPNSL);
    bool CALC_NSL = params.getBoolFlag(ARG_NSL);
    bool WRITE_DETAILED_IHS = params.getBoolFlag(ARG_IHS_DETAILED);
    bool CALC_XP = params.getBoolFlag(ARG_XP);
    bool CALC_SOFT = params.getBoolFlag(ARG_SOFT);
    bool SINGLE_EHH = false;

    this->MAF = params.getDoubleFlag(ARG_MAF);
    this->SKIP = params.getBoolFlag(ARG_SKIP);

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
                mapData.readMapDataVCF(vcfFilename, hapData.nloci);
            }
        }
        else
        {
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
        }
    }
    catch (...)
    {
        return 1;
    }

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

void MapData::readMapDataVCF(string filename, int expected_loci) {
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

    if (nloci != expected_loci)
    {
        cerr << "ERROR: Expected " << expected_loci << " loci in file but found " << nloci << ".\n";
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

    for (int i = 0; i < numCommentedLines; i++) {
        getline(fin, line);
    }

    this->initMapData(nloci); 
    
    double Mb = 1000000.0;
    
    string chr;
    for (int locus = 0; locus < this->nloci; locus++)
    {
        fin >> mapEntries[locus].chr;
        fin >> mapEntries[locus].physicalPos;
        fin >> mapEntries[locus].locusName;

        locus_query_map[mapEntries[locus].locusName] = locus;

        getline(fin, line);
        mapEntries[locus].geneticPos = double(mapEntries[locus].physicalPos)/Mb;
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
/*
void HapData::readHapData(string filename, bool unphased)
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
    int nhaps = 0;
    int previous_nloci = -1;
    int current_nloci = 0;
    //Counts number of haps (rows) and number of loci (cols)
    //if any lines differ, send an error message and throw an exception
    while (getline(fin, line))
    {
        //getline(fin,line);
        //if(fin.eof()) break;
        nhaps++;
        current_nloci = countFields(line);
        //cout << "nloci: " << current_nloci << endl;
        if (previous_nloci < 0)
        {
            previous_nloci = current_nloci;
            continue;
        }
        else if (previous_nloci != current_nloci)
        {
            cerr << "ERROR: line " << nhaps << " of " << filename << " has " << current_nloci
                 << ", but the previous line has " << previous_nloci << ".\n";
            throw 0;
        }
        previous_nloci = current_nloci;
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

    cerr << "Loading " << nhaps << " haplotypes and " << current_nloci << " loci...\n";
    

    if (unphased){
        initHapData(nhaps/2, current_nloci);
    }
    else{
        initHapData(nhaps, current_nloci);
    }

    char allele1;
    for (int hap = 0; hap < nhaps; hap++)
    {
        for (int locus = 0; locus < this->nloci; locus++)
        {
            if(unphased){
                fin >> allele1;
                if (allele1 != '0' && allele1 != '1'){
                    cerr << "ERROR: Alleles must be coded 0/1 only.\n";
                    cerr << allele1 << endl;
                    throw 0;
                }

                if (hap % 2 == 1){
                    if (allele1 == '1' && data->data[(hap-1)/2][locus] == '1'){
                        data->data[(hap-1)/2][locus] = '2';
                    }
                    else if (allele1 == '1' && data->data[(hap-1)/2][locus] == '0') {
                        data->data[(hap-1)/2][locus] = '1';
                    }
                }
                else{
                    data->data[hap/2][locus] = allele1;
                }
            }
            else{
                fin >> data->data[hap][locus];
                if (data->data[hap][locus] != '0' && data->data[hap][locus] != '1')
                {
                    cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
                    throw 0;
                }
            }
        }
    }

    fin.close();
}
*/

void HapData::readHapDataTPED(string filename, bool unphased)
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
    int nloci = 0;
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
                }
                else if (allele1 == '1' || allele2 == '1'){
                    allele = '1';
                    hapEntries[locus].positions.push_back(hap);
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
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    int numMapCols = 9;
    //int fileStart = fin.tellg();
    string line;
    int nloci = 0;
    int previous_nhaps = -1;
    int current_nhaps = 0;
    //Counts number of haps (cols) and number of loci (rows)
    //if any lines differ, send an error message and throw an exception
    while (getline(fin, line))
    {
        if (line[0] == '#') {
            continue;
        }
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

    int nhaps = unphased ? (current_nhaps - numMapCols) : (current_nhaps - numMapCols) * 2;
    int nfields = (current_nhaps - numMapCols);
    
    
    if(SKIP){
        cerr << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
        (*flog)  << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
    
    }

    cerr << "Loading " << nhaps << " haplotypes and " << nloci << " loci...\n";

    initHapData(nhaps, nloci);

    string junk;
    char allele1, allele2, separator;
    bool skipLine = false;
    for (int locus = 0; locus < this->nloci; locus++)
    {
        for (int i = 0; i < numMapCols; i++) {
            fin >> junk;
            if (i == 0 && junk[0] == '#') {
                skipLine = true;
                break;
            }
        }
        if (skipLine) {
            getline(fin, junk);
            skipLine = false;
            locus--;
            continue;
        }
        for (int field = 0; field < nfields; field++)
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
                    allele = '2';
                }
                else if (allele1 == '1' || allele2 == '1'){
                    allele = '1';
                }
                data->data[field][locus] = allele;
            }
            else{
                data->data[2 * field][locus] = allele1;
                data->data[2 * field + 1][locus] = allele2;
            }
        }
    }
 
    if(SKIP){
         cerr << "Removed " << mapData->nloci - count << " low frequency variants.\n";
        (*flog) << "Removed " << mapData->nloci - count << " low frequency variants.\n";
    }
   

    fin.close();
}


void HapData::initHapData(unsigned int nhaps, unsigned int nloci)
{
    if (nhaps < 1 || nloci < 1)
    {
        cerr << "ERROR: number of haplotypes (" << nhaps << ") and number of loci (" << nloci << ") must be positive.\n";
        throw 0;
    }

    this->hapEntries = new struct HapEntry[nhaps];

    
    this->nhaps = nhaps;
    this->nloci = nloci;

    //this->data = new char *[nhaps];
    // for (unsigned int i = 0; i < nhaps; i++)
    // {
    //     data->data[i] = new char[nloci];
    //     for (unsigned int j = 0; j < nloci; j++)
    //     {
    //         data->data[i][j] = MISSING_CHAR;
    //     }
    // }
    for (unsigned int j = 0; j < nloci; j++)
    {
        hapEntries[j].xors.reserve(nhaps);
        hapEntries[j].positions.reserve(nhaps);
    }
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

