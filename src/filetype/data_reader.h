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


#endif