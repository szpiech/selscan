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

#include <thread>

#include <ctime>
#include <cmath>

using namespace std;

class MainTools{
    protected:
        param_main p;
        
    public:
        HapMap hm;
        ofstream* flog;
        ofstream* fout;
        Bar *bar;
        int numThreads;
        MainTools(HapMap& hm, param_main& params,  ofstream* flog,  ofstream* fout);

        double static readTimer() {
            return clock() / (double) CLOCKS_PER_SEC;
        }
        inline unsigned int square_alt(int n){
            return n*n;
        }
        inline unsigned int twice_num_pair(int n){
            return n*n - n;
        }
        inline unsigned int num_pair(int n){
            return (n*n - n)/2;
        }

        int physicalDistance(int currentLocus, bool downstream){
            int distance;
            if(downstream){
                distance =  hm.mapData.mapEntries[currentLocus+1].physicalPos - hm.mapData.mapEntries[currentLocus].physicalPos;
            }else{
                distance = hm.mapData.mapEntries[currentLocus].physicalPos - hm.mapData.mapEntries[currentLocus-1].physicalPos;
            }

            // this should not happen as we already did integrity check previously
            if (distance < 0)
            {
                std::cerr << "ERROR: physical position not in ascending order.\n"; 
                throw 0;
            }
            return distance;
        }
        double geneticDistance(int currentLocus, bool downstream){
            double distance;
            if(downstream){
                distance =  hm.mapData.mapEntries[currentLocus+1].geneticPos - hm.mapData.mapEntries[currentLocus].geneticPos;
            }else{
                distance = hm.mapData.mapEntries[currentLocus].geneticPos - hm.mapData.mapEntries[currentLocus-1].geneticPos;
            }

            // this should not happen as we already did integrity check previously
            if (distance < 0)
            {
                std::cerr << "ERROR: genetic position not in ascending order.\n"; 
                throw 0;
            }
            return distance;
        }
    
    
};



class XPIHH : public MainTools{
    public:
        XPIHH(HapMap& hm, param_main& params,  ofstream* flog,  ofstream* fout) : MainTools(hm, params,  flog,  fout){        }
        void xpihh_main();
    
    private:
        double* ihh_p1;
        double* ihh_p2;
        void calc_xpihh(int locus);
        void calc_ehh_unidirection_xpihh(int locus, unordered_map<unsigned int, vector<unsigned int> > & m, bool downstream);
        void static thread_xpihh(int tid, unordered_map<unsigned int, vector<unsigned int> >& m, unordered_map<unsigned int, vector<unsigned int> >& md, XPIHH* obj);

};

class IHS : public MainTools{
    public:
        IHS(HapMap& hm, param_main& params,  ofstream* flog,  ofstream* fout) : MainTools(hm, params,  flog,  fout){        }
        void ihs_main(); //thread_ihs

    private:
        double* iHH0;
        double* iHH1;
        void static thread_ihs(int tid, unordered_map<unsigned int, vector<unsigned int> >& m, unordered_map<unsigned int, vector<unsigned int> >& md, IHS* ehh_obj);
        void calc_ehh_unidirection_ihs(int locus, unordered_map<unsigned int, vector<unsigned int> > & m, bool downstream);
        void calc_ihh(int locus);     
};

class EHH : public MainTools{
    public:
        EHH(HapMap& hm, param_main& params,  ofstream* flog,  ofstream* fout) : MainTools(hm, params,  flog,  fout){        }
        void calc_single_ehh(string query);
       
    private:
        void calc_ehh_unidirection(int locus, unordered_map<unsigned int, vector<unsigned int> > & m, bool downstream);
};



#endif