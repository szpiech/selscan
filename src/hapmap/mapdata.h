#ifndef __MAPDATA_H__
#define __MAPDATA_H__

#include <map>
#include <string>
#include <iostream>
#include <queue>

using namespace std;

struct MapEntry
{
    int physicalPos;
    double geneticPos;
    string locusName;
    string chr;
    int locId;
    bool skipLocus = false;
};

class MapData
{
public:
    struct MapEntry* mapEntries = NULL; //vector of map entries
    int nloci;
    map<string, int> locus_query_map;
    
    ~MapData(){
        releaseMapData();
    }    

    //allocates the arrays and populates them with -9 or "--" depending on type
    void initMapData(int nloci);
    void releaseMapData();

    //reads in map data and also does basic checks on integrity of format
    //returns a populated MapData structure if successful
    //throws an exception otherwise
    void readMapData(string filename, int expected_loci, bool USE_PMAP, queue<int>& skip_queue);
    void readMapDataTPED(string filename, int expected_loci, int expected_haps, bool USE_PMAP, queue<int>& skip_queue);
    void readMapDataVCF(string filename, int expected_loci, queue<int>& skipQueue); //Physical positions only

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

    void print(){
        for (int i = 0; i < nloci; i++)
        {
            cout << "Locus: " << i << " Physical Pos: " << mapEntries[i].physicalPos << endl;
            // cout << "Genetic Pos: " << mapEntries[i].geneticPos << endl;
            // cout << "Locus Name: " << mapEntries[i].locusName << endl;
            // cout << "Chromosome: " << mapEntries[i].chr << endl;
        }
    }

    /***
     * @param query: Locus name
     * @returns locus ( in range [0 to nloci) )
    */
    int queryFound(string query)
    {
        if (locus_query_map.count(query)>0){
            return locus_query_map[query];
        }
        return -1;
    }
};

#endif