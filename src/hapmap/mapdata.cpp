#include "mapdata.h"
#include <iostream>
#include <string>
#include "../gzstream.h"
#include <queue>
#include <sstream>

using namespace std;

/**
 * reads in map data and also does basic checks on integrity of format
 * @returns a populated MapData structure if successful
 * @throws an exception if not successful
 */
void MapData::readMapData(string filename, int expected_loci, bool USE_PMAP, queue<int>& skip_queue)
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
    int nloci_before_filter = 0;
    int num_cols = 4;
    int current_cols = 0;
    while (getline(fin, line))
    {
        nloci_before_filter++;
        current_cols = countFields(line);
        if (current_cols != num_cols)
        {
            cerr << "ERROR: line " << nloci_before_filter << " of " << filename << " has " << current_cols
                 << ", but expected " << num_cols << ".\n";
            throw 0;
        }
    }

    if (nloci_before_filter-skip_queue.size() != expected_loci)
    {
        cerr << "ERROR: Expected " << expected_loci << " loci in map file but found " << nloci_before_filter-skip_queue.size() << ".\n";
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

    cerr << "Loading map data for " << nloci_before_filter-skip_queue.size() << " loci\n";
    initMapData(nloci_before_filter-skip_queue.size());

    string chr;
    unsigned int locus_after_filter = 0;
    for (unsigned int locus_before_filter = 0; locus_before_filter < nloci_before_filter; locus_before_filter++)
    {
        if(!skip_queue.empty()){
            if(skip_queue.front()==locus_before_filter){
                skip_queue.pop();
                string junk;
                fin >> junk;
                fin >> junk;
                fin >>  junk;
                fin >>  junk;
                getline(fin, line);
                continue;
            }
        }

        fin >> mapEntries[locus_after_filter].chr;
        fin >> mapEntries[locus_after_filter].locusName;
        fin >> mapEntries[locus_after_filter].geneticPos;
        fin >> mapEntries[locus_after_filter].physicalPos;

        //locus_query_map[mapEntries[locus_after_filter].locusName] = locus_after_filter;
        mapEntries[locus_after_filter].locId = locus_before_filter;

        double Mb = 1000000.0;
        //if (USE_PMAP) 
        mapEntries[locus_after_filter].geneticPos = (long double)(mapEntries[locus_after_filter].physicalPos/Mb);
        //if (USE_PMAP) mapEntries[locus_after_filter].geneticPos = double(mapEntries[locus_after_filter].physicalPos);

        locus_after_filter++;
        getline(fin, line);
    }

    fin.close();
}

void MapData::readMapDataTPED(string filename, int expected_loci, int expected_haps, bool USE_PMAP, queue<int>& skip_queue)
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

    if (nloci - skip_queue.size() != expected_loci)
    {
        cerr << "ERROR: Expected " << expected_loci << " loci in map file but found " << nloci - skip_queue.size()  << ".\n";
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

    cerr << "Loading map data for " << nloci-skip_queue.size() << " loci\n";

    initMapData(nloci-skip_queue.size());
    
    double Mb = 1000000.0;
    
    string chr;
    int locus_after_filter = 0;
    for (int locus = 0; locus < this->nloci; locus++)  // locus = locus_before_filter
    {
        getline(fin, line);
        
        if(!skip_queue.empty()){
            if(skip_queue.front()==locus){
                skip_queue.pop();
                //getline(fin, line);
                continue;
            }
        }

        stringstream ss(line);
        ss >> mapEntries[locus_after_filter].chr;
        ss >> mapEntries[locus_after_filter].locusName;
        ss >> mapEntries[locus_after_filter].geneticPos;
        ss >> mapEntries[locus_after_filter].physicalPos;

        cout<<mapEntries[locus_after_filter].chr<<"\t" << mapEntries[locus_after_filter].locusName << "\t" << mapEntries[locus_after_filter].physicalPos << endl;
        //locus_query_map[mapEntries[locus_after_filter].locusName] = locus_after_filter;

        if (USE_PMAP) mapEntries[locus_after_filter].geneticPos = double(mapEntries[locus_after_filter].physicalPos)/Mb;
        //getline(fin, line);
        locus_after_filter++;
    }

    fin.close();
}

void MapData::readMapDataVCF(string filename, int expected_loci, queue <int>& skip_queue) {
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    string line;
    int nloci_before_filter = 0;
    int numCommentedLines = 0;
    while (getline(fin, line))
    {
        if (line[0] == '#') {
            numCommentedLines++;
        }
        else {
            nloci_before_filter++;
        }
    }

    if (nloci_before_filter-skip_queue.size() != expected_loci)
    {
        cerr << "ERROR: Expected " << expected_loci << " loci in file but found " << nloci_before_filter-skip_queue.size() << ".\n";
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

    cerr << "Loading map data for " << nloci_before_filter-skip_queue.size() << " loci\n";

    for (int i = 0; i < numCommentedLines; i++) {
        getline(fin, line);
    }

    this->initMapData(nloci_before_filter-skip_queue.size()); 
    
    double Mb = 1000000.0;
    
    string chr;

    int n_chromosomes_included = 0;
    unsigned int locus_after_filter = 0;
    for (unsigned  int locus_before_filter = 0; locus_before_filter < nloci_before_filter; locus_before_filter++)
    {
        if(!skip_queue.empty()){
            if(skip_queue.front()==locus_before_filter){
                skip_queue.pop();
                string junk;
                fin >> junk;
                fin >>  junk;
                fin >>  junk;
                getline(fin, line);
                continue;
            }
        }
        
        
        
        fin >> mapEntries[locus_after_filter].chr;
        fin >> mapEntries[locus_after_filter].physicalPos;
        fin >> mapEntries[locus_after_filter].locusName;
        //locus_query_map[mapEntries[locus].locusName] = locus;
        //locus_query_map[to_string(mapEntries[locus_after_filter].physicalPos)] = locus_after_filter;
        mapEntries[locus_after_filter].locId = locus_before_filter;


        //if exists in map do nothing else insert
        if(chr_list.find(mapEntries[locus_after_filter].chr) == chr_list.end()){
            chr_list[mapEntries[locus_after_filter].chr] = n_chromosomes_included++;
        }


        mapEntries[locus_after_filter].geneticPos = (long double)(mapEntries[locus_after_filter].physicalPos)/Mb;
        //cout<<mapEntries[locus_after_filter].geneticPos<<" "<<mapEntries[locus_after_filter].physicalPos<<endl;
        getline(fin, line);
        locus_after_filter++;
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

    mapEntries = new struct MapEntry[nloci];
    this->nloci = nloci;
}

void MapData::releaseMapData()
{
    if (mapEntries == NULL) return;
    this->nloci = -9;

    delete [] mapEntries;
    mapEntries = NULL;
    return;
}