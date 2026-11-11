#include "mapdata.h"
#include <iostream>
#include <string>
#include "../gzstream.h"
#include "../gzlinereader.h"
#include <queue>
#include <sstream>
#include <limits>

using namespace std;

/**
 * reads in map data and also does basic checks on integrity of format
 * @returns a populated MapData structure if successful
 * @throws an exception if not successful
 */
void MapData::readMapData(string filename, size_t expected_loci, bool USE_PMAP, queue<size_t>& skip_queue)
{
    igzstream fin;
    cerr << "Opening " << filename << " to read map data...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        exit(EXIT_FAILURE);
    }

    string line;
    size_t nloci_before_filter = 0;
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
            exit(EXIT_FAILURE);
        }
    }

    if (nloci_before_filter - (skip_queue.size()) != expected_loci)
    {
        cerr << "ERROR: Expected " << expected_loci << " loci in map file but found " << nloci_before_filter-skip_queue.size() << ".\n";
        exit(EXIT_FAILURE);
    }

    fin.clear(); // clear error flags
    fin.close();


    fin.open(filename.c_str());
    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        exit(EXIT_FAILURE);
    }

    cerr << "Loading map data for " << nloci_before_filter-skip_queue.size() << " loci\n";
    initMapData(nloci_before_filter-skip_queue.size());

    string chr;
    size_t locus_after_filter = 0;
    for (size_t locus_before_filter = 0; locus_before_filter < nloci_before_filter; locus_before_filter++)
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

       // mapEntries[locus_after_filter].geneticPos = (long double)(mapEntries[locus_after_filter].physicalPos/Mb);
        if (USE_PMAP) mapEntries[locus_after_filter].geneticPos = double(mapEntries[locus_after_filter].physicalPos/Mb);

        locus_after_filter++;
        getline(fin, line);
    }

    fin.close();
}

void MapData::readMapDataTPED(string filename, size_t expected_loci, int expected_haps, bool USE_PMAP, queue<size_t>& skip_queue)
{   
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        exit(EXIT_FAILURE);
    }

    string line;
    size_t nloci = 0;
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
            exit(EXIT_FAILURE);
        }
    }

    if (nloci - (skip_queue.size()) != expected_loci)
    {
        cerr << "ERROR: Expected " << expected_loci << " loci in map file but found " << nloci - skip_queue.size()  << ".\n";
        exit(EXIT_FAILURE);
    }

    fin.clear(); // clear error flags
    fin.close();



    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        exit(EXIT_FAILURE);
    }

    cerr << "Loading map data for " << nloci-skip_queue.size() << " loci\n";

    initMapData(nloci-skip_queue.size());
    
    double Mb = 1000000.0;
    
    string chr;
    size_t locus_after_filter = 0;
    for (size_t locus = 0; locus < nloci; locus++)  // locus = locus_before_filter
    {
        if(!skip_queue.empty()){
            if(skip_queue.front()==locus){
                skip_queue.pop();
                getline(fin, line);
                continue;
            }
        }
        fin >> mapEntries[locus_after_filter].chr;
        fin >> mapEntries[locus_after_filter].locusName;
        fin >> mapEntries[locus_after_filter].geneticPos;
        fin >> mapEntries[locus_after_filter].physicalPos;
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip rest of the line

        if (USE_PMAP) mapEntries[locus_after_filter].geneticPos = double(mapEntries[locus_after_filter].physicalPos)/Mb;
    
        locus_after_filter++;
    }

    fin.close();
}

void MapData::readMapDataVCF(string filename, size_t expected_loci, queue <size_t>& skip_queue) {
    // P1: read with the buffered GzLineReader and take CHROM, POS and ID
    // from each line directly (igzstream with getline/>> was the slowest
    // part of VCF loading). Empty lines are skipped, as the VCF readers do.
    GzLineReader fin;
    cerr << "Opening " << filename << " to read map data...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        exit(EXIT_FAILURE);
    }

    string line;
    size_t nloci_before_filter = 0;
    while (getline(fin, line))
    {
        if (line.empty() || line[0] == '#') continue;
        nloci_before_filter++;
    }

    if (nloci_before_filter-skip_queue.size() != expected_loci)
    {
        cerr << "ERROR: Expected " << expected_loci << " loci in file but found " << nloci_before_filter-skip_queue.size() << ".\n";
        exit(EXIT_FAILURE);
    }

    fin.close();
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        exit(EXIT_FAILURE);
    }

    cerr << "Loading map data for " << nloci_before_filter-skip_queue.size() << " loci\n";

    this->initMapData(nloci_before_filter-skip_queue.size());

    const double Mb = 1000000.0;

    // next whitespace-delimited field of line starting at p (as `>>` would)
    auto nextField = [](const string &l, size_t &p) -> string {
        while (p < l.size() && (l[p] == ' ' || l[p] == '\t')) ++p;
        size_t q = p;
        while (q < l.size() && l[q] != ' ' && l[q] != '\t') ++q;
        string f = l.substr(p, q - p);
        p = q;
        return f;
    };

    size_t locus_before_filter = 0;
    size_t locus_after_filter = 0;
    while (locus_before_filter < nloci_before_filter && getline(fin, line))
    {
        if (line.empty() || line[0] == '#') continue;
        if (!skip_queue.empty() && skip_queue.front() == locus_before_filter) {
            skip_queue.pop();
            locus_before_filter++;
            continue;
        }
        size_t p = 0;
        MapEntry &e = mapEntries[locus_after_filter];
        e.chr = nextField(line, p);
        e.physicalPos = atoi(nextField(line, p).c_str());
        e.locusName = nextField(line, p);
        e.locId = locus_before_filter;
        e.geneticPos = (long double)(e.physicalPos)/Mb;
        locus_after_filter++;
        locus_before_filter++;
    }

    fin.close();
}

void MapData::initMapData(int nloci)
{
    if (nloci < 1)
    {
        cerr << "ERROR: number of loci (" << nloci << ") must be positive.\n";
        exit(EXIT_FAILURE);
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