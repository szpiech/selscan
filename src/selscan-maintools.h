/* selscan -- a program to calculate EHH-based scans for positive selection in genomes
   Copyright (C) 2014-2024  Zachary A Szpiech
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

#ifndef __SELSCAN_MAINTOOLS_H__
#define __SELSCAN_MAINTOOLS_H__

#include <string>
#include <map>
#include <cstdio>
#include "binom.h"
#include "param_t.h"
#include "selscan-data.h"
#include "selscan-pbar.h"
#include "hamming_t.h"
#include "selscan-cli.h"

#include <ctime>
#include <cmath>


using namespace std;



class IHHFinder{
    public:
        HapMap hm;
        int numHaps;
        int numSnps;
        double* iHH0;
        double* iHH1;
        ofstream* flog;
        ofstream* fout;
        Bar *bar;
        
    IHHFinder(HapMap& hm, param_t& params,  ofstream* flog,  ofstream* fout){
        this->flog = flog;
        this->fout = fout; 
        this->hm = hm;
        this->numHaps = hm.hapData.nhaps;
        this->numSnps = hm.mapData.nloci;
        this->SCALE_PARAMETER = params.getIntFlag(ARG_GAP_SCALE);
        this->MAX_GAP = params.getIntFlag(ARG_MAX_GAP);
        this->EHH_CUTOFF = params.getDoubleFlag(ARG_CUTOFF);
        this->ALT = params.getBoolFlag(ARG_ALT);
        this->WAGH = params.getBoolFlag(ARG_WAGH);
        this->TRUNC = params.getBoolFlag(ARG_TRUNC);
        this->MAF = params.getDoubleFlag(ARG_MAF);
        this->numThreads = params.getIntFlag(ARG_THREAD);
        this->MAX_EXTEND = ( params.getIntFlag(ARG_MAX_EXTEND) <= 0 ) ? hm.mapData.mapEntries[numSnps-1].physicalPos - hm.mapData.mapEntries[0].physicalPos : params.getIntFlag(ARG_MAX_EXTEND);
        this->SKIP = params.getBoolFlag(ARG_SKIP);
        this->WRITE_DETAILED_IHS = params.getBoolFlag(ARG_IHS_DETAILED);
        this->unphased = params.getBoolFlag(ARG_UNPHASED);
        //double (*calc)(map<string, int> &, int, bool) = p->calc;
        this->CALC_XPNSL = params.getBoolFlag(ARG_XPNSL);
        //int MAX_EXTEND;
        if (!CALC_XPNSL){
           // MAX_EXTEND = ( ARG_MAX_EXTEND <= 0 ) ? physicalPos[nloci - 1] - physicalPos[0] : p->params->getIntFlag(ARG_MAX_EXTEND);
        }
        else{
            MAX_EXTEND = ( params.getIntFlag(ARG_MAX_EXTEND_NSL) <= 0 ) ? hm.mapData.mapEntries[numSnps-1].geneticPos - hm.mapData.mapEntries[0].geneticPos : params.getIntFlag(ARG_MAX_EXTEND_NSL);
        }
    }
    ~IHHFinder(){
        // delete[] iHH0;
        // delete[] iHH1;
    }

    void calcSingleEHH(string query);
    void calc_xpihh(int id);
    void calc_ihs();

    private:
        int SCALE_PARAMETER, MAX_GAP;
        double EHH_CUTOFF;
        bool ALT,WAGH,TRUNC;
        double MAF;
        int numThreads, MAX_EXTEND;
        bool SKIP, WRITE_DETAILED_IHS, unphased;
        bool CALC_XPNSL;
        
        //double (*calc)(map<string, int> &, int, bool) = p->calc;

        void calc_EHH_unidirection(int locus, unordered_map<unsigned int, vector<unsigned int> > & m, bool downstream);
        

    inline unsigned int twice_num_pair(int n){
        return n*n - n;
    }

    inline unsigned int num_pair(int n){
        return (n*n - n)/2;
    }

    double static readTimer() {
        return clock() / (double) CLOCKS_PER_SEC;
    }

    
};


// struct work_order_t
// {
//     int queryLoc;
//     int id;

//     string filename;

//     HaplotypeData *hapData;

//     HaplotypeData *hapData1;
//     HaplotypeData *hapData2;

//     MapData *mapData;

//     double (*calc)(map<string, int> &, int, bool);

//     double *ihs;
//     double *ihhDerivedLeft;
//     double *ihhDerivedRight;
//     double *ihhAncestralLeft;
//     double *ihhAncestralRight;
//     double *freq;

//     double *ihh1;
//     double *freq1;

//     double *ihh2;
//     double *freq2;

//     ofstream *flog;
//     ofstream *fout;
//     Bar *bar;

//     param_t *params;
// };

struct triplet_t
{
    double h1;
    double h12;
    double h2dh1;
};

void calculatePi(HapMap &hm, int winsize, string outFilename);

triplet_t calculateSoft(map<string, int> &count, int total);

void query_locus(void *work_order);
void query_locus_soft(void *order);

void calc_ihs(HapMap &hm);
void calc_nsl(HapMap &hm);
void calc_xpihh(HapMap &hm);
void calc_soft_ihs(HapMap &hm);

//double calcFreq(int locus, bool unphased);
int queryFound(MapData *mapData, string query);
void fillColors(int **hapColor, map<string, int> &hapCount,
                string *haplotypeList, int hapListLength,
                int currentLoc, int &currentColor, bool left);
bool familyDidSplit(const string &hapStr, const int hapCount,
                    int **hapColor, const int nhaps, const int colorIndex,
                    const int previousLoc, string &mostCommonHap);

double calculateHomozygosity_Wagh(map<string, int> &count, int total, int derivedCount);
double calculateHomozygosity(map<string, int> &count, int total, bool ALT);



#endif