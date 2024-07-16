#ifndef __SELSCAN_STATS_H__
#define __SELSCAN_STATS_H__

#include <time.h>
#include "../binom.h"
#include "../selscan-data.h"

#ifdef __APPLE__
#include <TargetConditionals.h>
#if TARGET_OS_MAC
    #define IS_MACOS
    #include "../../../../../../Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include/pthread/pthread.h"
#endif
#endif


/// @brief Class to hold common stuff required to compute all statistics for selscan
class SelscanStats {
    public:
        const std::unique_ptr<HapMap>& hm;
        param_main p;
        ofstream* flog;
        ofstream* fout;
        int numThreads;

        SelscanStats(const std::unique_ptr<HapMap>& hm, param_main& p,  ofstream* flog,  ofstream* fout): hm(hm), p(p){
            this->flog = flog;
            this->fout = fout; 
            this->numThreads = p.numThreads;
        }

        ~SelscanStats(){

        }

        //virtual void main() = 0;

        double static readTimer() {
            return clock() / (double) CLOCKS_PER_SEC;
        }
        inline unsigned int square_alt(int n){
            return n*n;
        }

        inline long double twice_num_pair_or_square(int n, bool alt){
            if(n < 2){
                return 0;
            }
            if(alt){
                return nCk(n, 2)+n;
            }
            return 2*nCk(n, 2);
            //return n*n - n;
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
                distance =  hm->mapData->mapEntries[core].physicalPos - hm->mapData->mapEntries[currentLocus].physicalPos;
            }else{
                distance = hm->mapData->mapEntries[currentLocus].physicalPos - hm->mapData->mapEntries[core].physicalPos;
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
                if(currentLocus+1> hm->hapData->nloci-1){
                    std::cerr << "ERROR: wrong locus"<<endl;
                    throw 0;
                }
                distance =   hm->mapData->mapEntries[currentLocus+1].physicalPos - hm->mapData->mapEntries[currentLocus].physicalPos;
            }else{
                if(currentLocus-1<0){
                    std::cerr << "ERROR: wrong locus"<<endl;
                    throw 0;
                }
                distance =  hm->mapData->mapEntries[currentLocus].physicalPos - hm->mapData->mapEntries[currentLocus-1].physicalPos;
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
                distance =  hm->mapData->mapEntries[currentLocus+1].geneticPos - hm->mapData->mapEntries[currentLocus].geneticPos;
            }else{
                distance = hm->mapData->mapEntries[currentLocus].geneticPos - hm->mapData->mapEntries[currentLocus-1].geneticPos;
            }

            // this should not happen as we already did integrity check previously
            if (distance < 0)
            {
                std::cerr << "ERROR: genetic position not in ascending order.\n"; 
                throw 0;
            }
            return distance;
        }

        // double calcFreq(int locus)
        // {
        //     double total = 0;
        //     double freq = 0;

        //     //assuming no missing data in hapmap structure
            
        //     if(this->hm->hapData->unphased){
        //         if(p.LOW_MEM){
        //             throw "Not implemented";
        //             //freq = hapEntries[locus].hapbitset->count_1s();
        //             freq = this->hm->hapData->hapEntries[locus].hapbitset->num_1s;
        //         }
        //         freq = this->hm->hapData->hapEntries[locus].count1 + this->hm->hapData->hapEntries[locus].count2*2;
        //         //freq = hapEntries[locus].positions.size() + hapEntries[locus].positions2.size()*2;
        //         total = this->hm->hapData->nhaps*2;
                
        //     }else{
        //         freq = this->hm->hapData->hapEntries[locus].positions.size();
        //         if(p.LOW_MEM){
        //             //freq = hapEntries[locus].hapbitset->count_1s();
        //             freq = this->hm->hapData->hapEntries[locus].hapbitset->num_1s;
        //         }
                
        //         total = this->hm->hapData->nhaps;
        //     }
        //     return freq/total;
        // }
};

#endif