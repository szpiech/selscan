#ifndef __SELSCAN_PI_H__
#define __SELSCAN_PI_H__

#include "selscan-stats.h"
#include "../thread_pool.h"

using namespace std;

class PI: public SelscanStats{
    public:
        ofstream flog_obj;
        ofstream fout_obj;
        PI(const std::unique_ptr<HapMap>&  hm, param_main& params) : SelscanStats(hm, params){  
            flog = &flog_obj;
            fout = &fout_obj;
            init_flog_fout("pi");    
               
            if(p.UNPHASED){
                cerr<<("ERROR: Unphased data not supported for PI calculation");
                exit(1);
            } 
            this->winsize = p.PI_WIN;    
            if(this->winsize <= 0 ){
                cerr<<("ERROR: Invalid window size for PI calculation");
                exit(1);
            }
            if(this->winsize > (hm->mapData->mapEntries[ hm->mapData->nloci-1].physicalPos)){
                //this->winsize = hm->mapData->mapEntries[ hm->mapData->nloci-1].physicalPos-1;
                (*flog)<< "Window size for PI calculation is larger than the physical length of the chromosome."<<this->winsize<<endl;

            }
        }
        void main(); 
        //pair<double, double> calc_pi(int start, int end);  
        
    private:
        int winsize;
        //static pthread_mutex_t mutex_log;
};



#endif
