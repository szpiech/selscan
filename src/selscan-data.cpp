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

//reads in map data and also does basic checks on integrity of format
//returns a populated MapData structure if successful
//throws an exception otherwise
MapData *readMapData(string filename, int expected_loci)
{
    ifstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    int fileStart = fin.tellg();
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
    fin.seekg(fileStart);

    cerr << "Loading map data for " << nloci << " loci\n";

    MapData *data = initMapData(nloci);

    string chr;
    for (int locus = 0; locus < data->nloci; locus++)
    {
        fin >> data->chr;
        fin >> data->locusName[locus];
        fin >> data->geneticPos[locus];
        fin >> data->physicalPos[locus];
    }

    fin.close();
    return data;
}

MapData *readMapDataTPED(string filename, int expected_loci, int expected_haps)
{
    ifstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    int fileStart = fin.tellg();
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
    fin.seekg(fileStart);

    cerr << "Loading map data for " << nloci << " loci\n";

    MapData *data = initMapData(nloci);

    string chr;
    for (int locus = 0; locus < data->nloci; locus++)
    {
        fin >> data->chr;
        fin >> data->locusName[locus];
        fin >> data->geneticPos[locus];
        fin >> data->physicalPos[locus];
        getline(fin, line);
    }

    fin.close();
    return data;
}

//allocates the arrays and populates them with MISSING or "--" depending on type
MapData *initMapData(int nloci)
{
    if (nloci < 1)
    {
        cerr << "ERROR: number of loci (" << nloci << ") must be positive.\n";
        throw 0;
    }

    MapData *data = new MapData;
    data->nloci = nloci;
    data->locusName = new string[nloci];
    data->physicalPos = new int[nloci];
    data->geneticPos = new double[nloci];

    for (int locus = 0; locus < nloci; locus++)
    {
        data->locusName[locus] = "--";
        data->physicalPos[locus] = MISSING;
        data->geneticPos[locus] = MISSING;
    }

    return data;
}

void releaseMapData(MapData *data)
{
    if (data == NULL) return;
    data->nloci = -9;
    delete [] data->locusName;
    delete [] data->physicalPos;
    delete [] data->geneticPos;
    delete data;
    data = NULL;
    return;
}

//reads in haplotype data and also does basic checks on integrity of format
//returns a populated HaplotypeData structure if successful
//throws an exception otherwise
HaplotypeData *readHaplotypeData(string filename)
{
    ifstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    int fileStart = fin.tellg();
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
    fin.seekg(fileStart);

    cerr << "Loading " << nhaps << " haplotypes and " << current_nloci << " loci...\n";

    HaplotypeData *data = initHaplotypeData(nhaps, current_nloci);


    for (int hap = 0; hap < data->nhaps; hap++)
    {
        for (int locus = 0; locus < data->nloci; locus++)
        {
            fin >> data->data[hap][locus];
            if (data->data[hap][locus] != '0' && data->data[hap][locus] != '1')
            {
                cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
                throw 0;
            }
        }
    }

    fin.close();

    return data;
}

HaplotypeData *readHaplotypeDataTPED(string filename)
{
    ifstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    int numMapCols = 4;
    int fileStart = fin.tellg();
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
    fin.seekg(fileStart);

    cerr << "Loading " << current_nhaps - numMapCols << " haplotypes and " << nloci << " loci...\n";

    HaplotypeData *data = initHaplotypeData(current_nhaps - numMapCols, nloci);

    string junk;

    for (int locus = 0; locus < data->nloci; locus++)
    {
        for (int i = 0; i < numMapCols; i++)
        {
            fin >> junk;
        }
        for (int hap = 0; hap < data->nhaps; hap++)
        {
            fin >> data->data[hap][locus];
            if (data->data[hap][locus] != '0' && data->data[hap][locus] != '1')
            {
                cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
                throw 0;
            }
        }
    }

    fin.close();

    return data;
}

HaplotypeData *initHaplotypeData(unsigned int nhaps, unsigned int nloci)
{
    if (nhaps < 1 || nloci < 1)
    {
        cerr << "ERROR: number of haplotypes (" << nhaps << ") and number of loci (" << nloci << ") must be positive.\n";
        throw 0;
    }

    HaplotypeData *data = new HaplotypeData;
    data->nhaps = nhaps;
    data->nloci = nloci;

    data->data = new char *[nhaps];
    for (unsigned int i = 0; i < nhaps; i++)
    {
        data->data[i] = new char[nloci];
        for (unsigned int j = 0; j < nloci; j++)
        {
            data->data[i][j] = MISSING_CHAR;
        }
    }

    return data;
}

void releaseHapData(HaplotypeData *data)
{
    if (data == NULL) return;
    for (int i = 0; i < data->nhaps; i++)
    {
        delete [] data->data[i];
    }

    delete [] data->data;

    data->data = NULL;
    data->nhaps = -9;
    data->nloci = -9;
    data = NULL;
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

