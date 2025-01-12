#ifndef __SELSCAN_STATS_H__
#define __SELSCAN_STATS_H__

#include <time.h>
#include <thread>
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

    string get_filename_base(string statname){
        string outFilename = p.outFilename + "." + statname;
        if(statname == "pi"){
            outFilename += "." + to_string(p.PI_WIN) + "bp";
            //else if (p.CALC_PI) p.outFilename += ".pi." + string(PI_WIN_str) + "bp";
        }
        if(statname == "ehh" && p.SINGLE_EHH) outFilename += "." + p.query;
        if (p.ALT) outFilename += ".alt";
        return outFilename;
    }

    void init_multi(){
        //open multiple filestream for all item in chrlist
        ofstream* fouts[hm->mapData->chr_list.size()];
        for (const auto& [chr, i] : hm->mapData->chr_list){
            cout<<chr<<" " <<i<<endl;
            string filename = get_filename_base("ihs") + "." + chr ;
            fouts[i] = new ofstream(filename);
        }
    }

    void init_flog_fout(string statname){
        string outFilename = get_filename_base(statname);
        (*flog) << "\nv" + VERSION + "\nCalculating ";
        if (statname == "ihs") (*flog) << " iHS.\n";
        else if (statname == "nsl") (*flog) << " nSL.\n";
        else if (statname == "xpnsl") (*flog) << " XP-nSL.\n";
        else if (statname == "xpehh") (*flog) << " XP-EHH.\n";
        else if (statname == "ehh") (*flog) << " EHH.\n";
        else if (statname == "pi") (*flog) << " PI.\n";
        else if (statname == "ihh12") (*flog) << " iHH12.\n";

        if (p.TPED)
        {
            (*flog) << "Input filename: " << p.tpedFilename << "\n";
            if (p.CALC_XP || p.CALC_XPNSL) (*flog) << "Reference input filename: " << p.tpedFilename2 << "\n";
        }
        else if (p.VCF) {
            (*flog) << "Input filename: " << p.vcfFilename << "\n";
            if (p.CALC_XP || p.CALC_XPNSL) (*flog) << "Reference input filename: " << p.vcfFilename2 << "\n";
            (*flog) << "Map filename: " << p.mapFilename << "\n";
        }
        else {
            (*flog) << "Input filename: " << p.hapFilename << "\n";
            if (p.CALC_XP || p.CALC_XPNSL) (*flog) << "Reference input filename: " << p.hapFilename2 << "\n";
            (*flog) << "Map filename: " << p.mapFilename << "\n";
        }

        (*flog) << "Output file: " << p.outFilename << "\n";
        (*flog) << "Threads: " << p.numThreads << "\n";
        (*flog) << "Scale parameter: " << p.SCALE_PARAMETER << "\n";
        (*flog) << "Max gap parameter: " << p.MAX_GAP << "\n";
        (*flog) << "EHH cutoff value: " << p.EHH_CUTOFF << "\n";
                (*flog) << "Minimum maf cutoff: " << p.MAF << "\n";
        (*flog) << "Missing value: " << p.MISSING_ALLOWED << "\n";

        (*flog) << "Phased: ";
        if(p.UNPHASED) (*flog) << "no\n";
        else (*flog) << "yes\n";
        (*flog) << "Alt flag set: ";
        if (p.ALT) (*flog) << "yes\n";
        else (*flog) << "no\n";


        if (p.ALT) (*flog) << "Missing handle mode: \n";
        else (*flog) << p.MISSING_MODE << "\n";
        (*flog).flush();
        
    }
        //SelscanStats(const std::unique_ptr<HapMap>& hm, param_main& p,  ofstream* flog,  ofstream* fout): hm(hm), p(p){
        SelscanStats(const std::unique_ptr<HapMap>& hm, param_main& p): hm(hm), p(p){
            this->numThreads = p.numThreads;
            string outFilename = p.outFilename;
            this->flog = hm->flog;
        }

        ~SelscanStats(){
            //flog->close();
        }
        

        double static readTimer() {
            return clock() / (double) CLOCKS_PER_SEC;
        }
        
        inline double square_alt(int n){
            return n*n;
        }

        inline double twice_num_pair_or_square(int n, bool alt){
            if(n < 2){
                return 0;
            }
            if(alt){
                return nCk(n, 2)+n;
            }
            return 2*nCk(n, 2);
        }

        inline long double twice_num_pair(int n){
            if(n < 2){
                return 0;
            }
            return 2*nCk(n, 2);
        }

        inline double num_pair(int n){
            if(n < 2){
                return 0;
            }
            return nCk(n, 2);
        }


        double geneticDistance_from_core(int currentLocus, int core, bool downstream){
            double distance;
            if(downstream){
                distance =  hm->mapData->mapEntries[core].geneticPos - hm->mapData->mapEntries[currentLocus].geneticPos;
            }else{
                distance = hm->mapData->mapEntries[currentLocus].geneticPos - hm->mapData->mapEntries[core].geneticPos;
            }

            // this should not happen as we already did integrity check previously
            if (distance < 0)
            {
                cout<<"Distance: current: core: isLeft:"<<distance<<" "<<currentLocus<<" "<<core<<" "<<downstream<<"\n";
                std::cerr << "ERROR: physical position not in ascending order.\n"; 
                throw 0;
            }
            return distance;
        }

        double physicalDistance_from_core(int currentLocus, int core, bool downstream){
            int distance;
            if(downstream){
                distance =  hm->mapData->mapEntries[core].physicalPos - hm->mapData->mapEntries[currentLocus].physicalPos;
            }else{
                distance = hm->mapData->mapEntries[currentLocus].physicalPos - hm->mapData->mapEntries[core].physicalPos;
            }

            // this should not happen as we already did integrity check previously
            if (distance < 0)
            {
                cout<<"Distance: current: core: isLeft:"<<distance<<" "<<currentLocus<<" "<<core<<" "<<downstream<<"\n";
                std::cerr << "ERROR: physical position not in ascending order.\n"; 
                throw 0;
            }
            return distance;
        }

        /// @brief Compute physical distance between currentLocus and previous one
        double physicalDistance(int currentLocus, bool downstream){
            double distance;
            if(downstream){
                if(currentLocus+1> hm->hapData->nloci-1){
                    std::cerr << "ERROR: wrong locus"<<endl;
                    throw 0;
                }
                distance =   double(hm->mapData->mapEntries[currentLocus+1].physicalPos - hm->mapData->mapEntries[currentLocus].physicalPos);
            }else{
                if(currentLocus-1<0){
                    std::cerr << "ERROR: wrong locus"<<endl;
                    throw 0;
                }
                distance =  double(hm->mapData->mapEntries[currentLocus].physicalPos - hm->mapData->mapEntries[currentLocus-1].physicalPos);
                // if(currentLocus == 22){
                //     cout<<hm->mapData->mapEntries[currentLocus].physicalPos<<" "<<hm->mapData->mapEntries[currentLocus-1].physicalPos<<" "<<distance<<endl;
                // }
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
            
        /// @brief Compute genetic distance between currentLocus and previous one    
        double geneticDistance(int currentLocus, bool downstream){
            double distance;
            if(downstream){
                distance =  double(hm->mapData->mapEntries[currentLocus+1].geneticPos - hm->mapData->mapEntries[currentLocus].geneticPos);
            }else{
                distance = double(hm->mapData->mapEntries[currentLocus].geneticPos - hm->mapData->mapEntries[currentLocus-1].geneticPos);
            }

            // this should not happen as we already did integrity check previously
            if (distance < 0)
            {
                std::cerr << "ERROR: genetic position not in ascending order.\n"; 
                throw 0;
            }
            return distance;
        }


        /// @brief Display progress bar
        void static displayProgressBar(std::atomic<int>& current, int total) {
            int barWidth = 70;  // Width of the progress bar
            while (current < total) {
                std::cout << "[";
                int pos = barWidth * current / total;
                for (int i = 0; i < barWidth; ++i) {
                    if (i < pos) std::cout << "=";
                    else if (i == pos) std::cout << ">";
                    else std::cout << " ";
                }
                std::cout << "] " << int(current * 100.0 / total) << " % (" << current << "/" << total << ")\r";
                std::cout.flush();
                std::this_thread::sleep_for(std::chrono::milliseconds(100));  // Update interval
            }
            // Final display to complete the progress bar at 100%
            std::cout << "[";
            for (int i = 0; i < barWidth; ++i) {
                std::cout << "=";
            }
            std::cout << "] 100 % (" << total << "/" << total << ")" << std::endl;
        }

        bool nextLocOutOfBounds(int locus, bool downstream){
            if(!downstream){
                if(locus+1 >= hm->hapData->nloci){
                    return true;
                }
                if(hm->mapData->mapEntries[locus+1].chr != hm->mapData->mapEntries[locus].chr){
                    cerr<<"Chr-out-of-bounds: pos"<<locus<<endl;
                    return true;
                }
            }else{
                if(locus-1 < 0){
                    return true;
                }
                if(hm->mapData->mapEntries[locus-1].chr != hm->mapData->mapEntries[locus].chr){
                    cerr<<"Chr-out-of-bounds: pos"<<locus<<endl;
                    return true;
                }
            }
            return false;
        }

        int getChrIdxFromLoc(int locus){
            return hm->mapData->chr_list[hm->mapData->mapEntries[locus].chr];
        }
};

#endif