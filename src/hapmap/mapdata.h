#ifndef __MAPDATA_H__
#define __MAPDATA_H__

#include <map>
#include <string>
#include <iostream>
#include <queue>

using namespace std;

// B44: physical positions are stored as int. A position field must be a
// whole number from 0 to 2,147,483,647; anything else (including the
// positions of chromosomes longer than 2.1 Gb) is rejected. 3.0.1 wrapped
// or clamped such values silently.
const char* const POS_LIMIT_MSG = " is not a whole number from 0 to 2147483647. Positions of 2^31 bp or more are not supported in this version.";
inline bool parsePosition(const char* s, const char* e, int &out)
{
    if (s >= e) return false;
    long long x = 0;
    for (const char* p = s; p < e; ++p) {
        if (*p < '0' || *p > '9') return false;
        x = x * 10 + (*p - '0');
        if (x > 2147483647LL) return false;
    }
    out = static_cast<int>(x);
    return true;
}
inline bool parsePosition(const string &f, int &out) { return parsePosition(f.data(), f.data() + f.size(), out); }

struct MapEntry
{
    int physicalPos; // uint64_t physicalPos;
    double geneticPos;
    string locusName;
    string chr = ".";
    int locId;
    //bool skipLocus = false;
};

class MapData
{
public:
    struct MapEntry* mapEntries = NULL; //vector of map entries
    int nloci;
    
    //map<string, int> chr_list;
    
    //map<string, int> locus_query_map;
    
    ~MapData(){
        releaseMapData();
    }    

    //allocates the arrays and populates them with -9 or "--" depending on type
    void initMapData(int nloci);
    void setFromEntries(const vector<MapEntry> &entries);
    void releaseMapData();

    //reads in map data and also does basic checks on integrity of format
    //returns a populated MapData structure if successful
    //throws an exception otherwise
    void readMapData(string filename, size_t expected_loci, bool USE_PMAP, queue<size_t>& skip_queue);
    void readMapDataTPED(string filename, size_t expected_loci, int expected_haps, bool USE_PMAP, queue<size_t>& skip_queue);
    void readMapDataVCF(string filename, size_t expected_loci, queue<size_t>& skipQueue); //Physical positions only

    inline int countFields(const string &str)
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

    /// FOR_DEBUGGING
    void print(int loc=5){
        //loc = nloci
        for (int i = 0; i < loc; i++)
        {
            //cout<<i <<" "<<mapEntries[i].physicalPos<<" "<<mapEntries[i].geneticPos<<" "<<mapEntries[i].locusName<<endl;
            cout << "Locus: " << i << " Physical Pos: " << mapEntries[i].physicalPos << " LocId: " << mapEntries[i].locId << endl;
            // cout << "Genetic Pos: " << mapEntries[i].geneticPos << endl;
            // cout << "Locus Name: " << mapEntries[i].locusName << endl;
            // cout << "Chromosome: " << mapEntries[i].chr << endl;
        }
    }
};

#endif