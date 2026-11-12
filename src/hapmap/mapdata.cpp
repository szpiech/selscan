#include "mapdata.h"
#include <iostream>
#include <string>
#include "../gzstream.h"
#include "../gzlinereader.h"
#include <queue>
#include <sstream>
#include <limits>
#include <cctype>

using namespace std;

// true if the line has a non-whitespace character (U9: blank lines are skipped)
static inline bool hasData(const string &l)
{
    for (char c : l) if (!std::isspace(static_cast<unsigned char>(c))) return true;
    return false;
}

// U1: the map and the genotype input must have the same number of loci
// before filtering. skip_queue holds the genotype loci removed by filters,
// so expected_loci + skip_queue.size() is the genotype file's count.
static void checkLocusCount(const string &filename, size_t mapRows, size_t expected_loci, size_t nskipped, const string &what)
{
    if (mapRows != expected_loci + nskipped) {
        cerr << "ERROR: " << filename << " has " << mapRows << " " << what << ", but the genotype data have "
             << expected_loci + nskipped << " loci (before filtering). They must match one to one.\n";
        exit(EXIT_FAILURE);
    }
}

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
    size_t lineNum = 0;
    while (getline(fin, line))
    {
        lineNum++;
        if (!hasData(line)) continue;   // U9
        nloci_before_filter++;
        current_cols = countFields(line);
        if (current_cols != num_cols)
        {
            cerr << "ERROR: line " << lineNum << " of " << filename << " has " << current_cols
                 << " columns, but expected " << num_cols << ".\n";
            exit(EXIT_FAILURE);
        }        {   // B44: column 4 is the physical position
            istringstream ss(line); string a, b, c, f4; int d;
            ss >> a >> b >> c >> f4;
            if (!parsePosition(f4, d)) {
                cerr << "ERROR: position '" << f4 << "' on line " << lineNum << " of " << filename << POS_LIMIT_MSG << "\n";
                exit(EXIT_FAILURE);
            }
        }
    }

    checkLocusCount(filename, nloci_before_filter, expected_loci, skip_queue.size(), "rows");

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
        while (getline(fin, line) && !hasData(line)) {}   // U9
        if(!skip_queue.empty()){
            if(skip_queue.front()==locus_before_filter){
                skip_queue.pop();
                continue;
            }
        }

        istringstream ss(line);
        ss >> mapEntries[locus_after_filter].chr;
        ss >> mapEntries[locus_after_filter].locusName;
        ss >> mapEntries[locus_after_filter].geneticPos;
        ss >> mapEntries[locus_after_filter].physicalPos;
        

        //locus_query_map[mapEntries[locus_after_filter].locusName] = locus_after_filter;
        mapEntries[locus_after_filter].locId = locus_before_filter;

        double Mb = 1000000.0;

       // mapEntries[locus_after_filter].geneticPos = (long double)(mapEntries[locus_after_filter].physicalPos/Mb);
        if (USE_PMAP) mapEntries[locus_after_filter].geneticPos = double(mapEntries[locus_after_filter].physicalPos/Mb);

        locus_after_filter++;
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
    size_t lineNum = 0;
    while (getline(fin, line))
    {
        lineNum++;
        if (!hasData(line)) continue;   // U9
        nloci++;
        current_cols = countFields(line);
        if (current_cols != num_cols + expected_haps)
        {
            cerr << "ERROR: line " << lineNum << " of " << filename << " has " << current_cols
                 << " columns, but expected " << num_cols + expected_haps << ".\n";
            exit(EXIT_FAILURE);
        }
    }

    checkLocusCount(filename, nloci, expected_loci, skip_queue.size(), "rows");

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
        while (getline(fin, line) && !hasData(line)) {}   // U9
        if(!skip_queue.empty()){
            if(skip_queue.front()==locus){
                skip_queue.pop();
                continue;
            }
        }
        istringstream ss(line);
        ss >> mapEntries[locus_after_filter].chr;
        ss >> mapEntries[locus_after_filter].locusName;
        ss >> mapEntries[locus_after_filter].geneticPos;
        ss >> mapEntries[locus_after_filter].physicalPos;

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

    checkLocusCount(filename, nloci_before_filter, expected_loci, skip_queue.size(), "variant records");

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

// P1: install map entries collected while the VCF haplotypes were read
void MapData::setFromEntries(std::deque<MapEntry> &entries)
{
    const size_t n = entries.size();
    cerr << "Loading map data for " << n << " loci (positions from the VCF)\n";
    this->initMapData((int)n);
    for (size_t i = 0; i < n; i++) {
        mapEntries[i] = std::move(entries.front());
        entries.pop_front();
    }
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