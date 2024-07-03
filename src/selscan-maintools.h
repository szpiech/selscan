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


// #include "stats/ihs.h"

# include <unordered_map>
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
        inline double twice_num_pair(int n){
            if(n < 2){
                return 0;
            }
            return 2*nCk(n, 2);
            //return n*n - n;
        }
        inline unsigned int num_pair(int n){
            return (n*n - n)/2;
        }



        int physicalDistance_from_core(int currentLocus, int core, bool downstream){
            int distance;
            if(downstream){
                distance =  hm.mapData.mapEntries[core].physicalPos - hm.mapData.mapEntries[currentLocus].physicalPos;
            }else{
                distance = hm.mapData.mapEntries[currentLocus].physicalPos - hm.mapData.mapEntries[core].physicalPos;
            }

            // this should not happen as we already did integrity check previously
            if (distance < 0)
            {
                cout<<"Distance: "<<distance<<" "<<currentLocus<<" "<<downstream<<"\n";
                std::cerr << "ERROR: physical position not in ascending order.\n"; 
                throw 0;
            }
            return distance;
        }

        int physicalDistance(int currentLocus, bool downstream){
            int distance;
            if(downstream){
                if(currentLocus+1>hm.hapData.nloci-1){
                    std::cerr << "ERROR: wrong locus"<<endl;
                    throw 0;
                }
                distance =  hm.mapData.mapEntries[currentLocus+1].physicalPos - hm.mapData.mapEntries[currentLocus].physicalPos;
            }else{
                if(currentLocus-1<0){
                    std::cerr << "ERROR: wrong locus"<<endl;
                    throw 0;
                }
                distance = hm.mapData.mapEntries[currentLocus].physicalPos - hm.mapData.mapEntries[currentLocus-1].physicalPos;
            }

            // this should not happen as we already did integrity check previously
            if (distance < 0)
            {
                cout<<"Distance: "<<distance<<" "<<currentLocus<<" "<<downstream<<"\n";
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





#endif