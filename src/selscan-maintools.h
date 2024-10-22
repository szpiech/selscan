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
#include "selscan-cli.h"

#include "stats/ihs.h"
#include "stats/ehh.h"
#include "stats/xpihh.h"
#include "stats/ihh12.h"
#include "stats/pi.h"


#include <unordered_map>
#include <thread>

#include <ctime>
#include <cmath>

using namespace std;

class MainTools{
    protected:
        param_main p;
        std::unique_ptr<HapMap> hm; 
        
    public:
        MainTools(param_t& params){
            p.hapFilename = params.getStringFlag(ARG_FILENAME_POP1);
            p.hapFilename2 = params.getStringFlag(ARG_FILENAME_POP2);
            p.mapFilename = params.getStringFlag(ARG_FILENAME_MAP);
            p.tpedFilename = params.getStringFlag(ARG_FILENAME_POP1_TPED);
            p.tpedFilename2 = params.getStringFlag(ARG_FILENAME_POP2_TPED);
            p.vcfFilename = params.getStringFlag(ARG_FILENAME_POP1_VCF);
            p.vcfFilename2 = params.getStringFlag(ARG_FILENAME_POP2_VCF);
            
            p.TPED = false;
            if (p.tpedFilename.compare(DEFAULT_FILENAME_POP1_TPED) != 0) p.TPED = true;

            p.VCF = false;
            if (p.vcfFilename.compare(DEFAULT_FILENAME_POP1_VCF) != 0) p.VCF = true;

            p.outFilename = params.getStringFlag(ARG_OUTFILE);
            p.query = params.getStringFlag(ARG_EHH);

            p.numThreads = params.getIntFlag(ARG_THREAD);
            p.SCALE_PARAMETER = params.getIntFlag(ARG_GAP_SCALE);
            p.MAX_GAP = params.getIntFlag(ARG_MAX_GAP);

            p.EHH_CUTOFF = params.getDoubleFlag(ARG_CUTOFF);
            p.MAF = params.getDoubleFlag(ARG_MAF);

            p.UNPHASED = params.getBoolFlag(ARG_UNPHASED);
            p.USE_PMAP = params.getBoolFlag(ARG_PMAP);
            p.ALT = params.getBoolFlag(ARG_ALT);
            p.WAGH = params.getBoolFlag(ARG_WAGH);
            p.CALC_IHS = params.getBoolFlag(ARG_IHS);
            p.CALC_XPNSL = params.getBoolFlag(ARG_XPNSL);
            p.CALC_NSL = params.getBoolFlag(ARG_NSL);
            p.WRITE_DETAILED_IHS = params.getBoolFlag(ARG_IHS_DETAILED);
            p.CALC_XP = params.getBoolFlag(ARG_XP);
            p.CALC_SOFT = params.getBoolFlag(ARG_SOFT);
            p.SINGLE_EHH = false;

            p.SKIP = !params.getBoolFlag(ARG_KEEP);//params.getBoolFlag(ARG_SKIP);
            if(params.getBoolFlag(ARG_SKIP)){
                cerr << "WARNING: " << ARG_SKIP << " is now on by dafault.  This flag no longer has a function.\n";
            }
            
            // p.EHH1K = params.getIntFlag(ARG_SOFT_K);
            p.QWIN = params.getIntFlag(ARG_QWIN);

            p.CALC_PI = params.getBoolFlag(ARG_PI);
            p.PI_WIN = params.getIntFlag(ARG_PI_WIN);

            p.LOW_MEM = params.getBoolFlag(ARG_LOW_MEM);
            p.MISSING = params.getBoolFlag(ARG_MISSING_FLAG);

            
            p.MAX_EXTEND_NSL = params.getIntFlag(ARG_MAX_EXTEND_NSL); //p.MAX_EXTEND_NSL = ( params.getIntFlag(ARG_MAX_EXTEND_NSL) <= 0 ) ? hm.mapData.nloci : params.getIntFlag(ARG_MAX_EXTEND_NSL);
            p.MAX_EXTEND =  params.getIntFlag(ARG_MAX_EXTEND);
            p.TRUNC = params.getBoolFlag(ARG_TRUNC);

            
            // this->hm = hm;
            //this->bar = new Bar();
            //this->numThreads = p.numThreads
        }
        
        ~MainTools(){
            //delete hm;
        }

        int run(int argc, char *argv[]){  
            auto start = std::chrono::high_resolution_clock::now();
            double time_start = (clock() / (double) CLOCKS_PER_SEC);

            char PI_WIN_str[50];
            snprintf(PI_WIN_str,50, "%d", p.PI_WIN);

            if (p.query.compare(DEFAULT_EHH) != 0) p.SINGLE_EHH = true;

            if (p.SINGLE_EHH) p.outFilename += ".ehh." + p.query;
            else if (p.CALC_IHS) p.outFilename += ".ihs";
            else if (p.CALC_NSL) p.outFilename += ".nsl";
            else if (p.CALC_XPNSL) p.outFilename += ".xpnsl";
            else if (p.CALC_XP) p.outFilename += ".xpehh";
            else if (p.CALC_SOFT) p.outFilename += ".ihh12";
            else if (p.CALC_PI) p.outFilename += ".pi." + string(PI_WIN_str) + "bp";

            if (p.ALT) p.outFilename += ".alt";

            //Open stream for log file
            ofstream flog;
            string logFilename = p.outFilename + ".log";
            flog.open(logFilename.c_str());
            if (flog.fail())
            {
                cerr << "ERROR: could not open " << logFilename << " for writing.\n";
                return 1;
            }

            //Open stream for output file
            ofstream fout;
            p.outFilename += ".out";
            fout.open(p.outFilename.c_str());
            if (fout.fail())
            {
                cerr << "ERROR: could not open " << p.outFilename << " for writing.\n";
                return 1;
            }

            // input data is loaded into HapMap object
            //hm = new HapMap(p, &flog, &fout);
            //std::unique_ptr<HapMap> 

            hm = std::make_unique<HapMap>(p, &flog, &fout);
            bool ERROR = hm->loadHapMapData();
            if (ERROR)
            {
                cerr << "ERROR: could not load hapmap data.\n";
                throw 1;
            }

            for (int i = 0; i < argc; i++)
            {
                flog << argv[i] << " ";
            }

            flog << "\nv" + VERSION + "\nCalculating ";
            if (p.CALC_XP) flog << "XP-EHH.\n";
            else if (p.CALC_PI) flog << "PI.\n";
            else if (p.CALC_IHS) flog << " iHS.\n";
            else if (p.CALC_NSL) flog << " nSL.\n";
            else if (p.CALC_XPNSL) flog << " XP-nSL.\n";
            else if (p.CALC_SOFT) flog << " iHH1K.\n";

            if (p.TPED)
            {
                flog << "Input filename: " << p.tpedFilename << "\n";
                if (p.CALC_XP || p.CALC_XPNSL) flog << "Reference input filename: " << p.tpedFilename2 << "\n";
            }
            else if (p.VCF) {
                flog << "Input filename: " << p.vcfFilename << "\n";
                if (p.CALC_XP || p.CALC_XPNSL) flog << "Reference input filename: " << p.vcfFilename2 << "\n";
                flog << "Map filename: " << p.mapFilename << "\n";
            }
            else {
                flog << "Input filename: " << p.hapFilename << "\n";
                if (p.CALC_XP || p.CALC_XPNSL) flog << "Reference input filename: " << p.hapFilename2 << "\n";
                flog << "Map filename: " << p.mapFilename << "\n";
            }
            flog << "Output file: " << p.outFilename << "\n";
            flog << "Threads: " << p.numThreads << "\n";
            flog << "Scale parameter: " << p.SCALE_PARAMETER << "\n";
            flog << "Max gap parameter: " << p.MAX_GAP << "\n";
            flog << "EHH cutoff value: " << p.EHH_CUTOFF << "\n";
            flog << "Missing value: " << p.MISSING << "\n";

            flog << "Phased: ";
            if(p.UNPHASED) flog << "no\n";
            else flog << "yes\n";
            flog << "Alt flag set: ";
            if (p.ALT) flog << "yes\n";
            else flog << "no\n";
            flog.flush();

            if (hm->mapData->nloci < p.numThreads)
            {
                p.numThreads = 1;
                cerr << "WARNING: there are fewer loci than threads requested.  Running with " << p.numThreads << " thread instead.\n";
            }

            if (p.SINGLE_EHH)
            {
                EHH ehhfinder(hm, p, &flog, &fout);
                if (p.CALC_SOFT)
                {
                    throw ("ERROR: Soft EHH not implemented yet.\n");
                }
                else
                {
                    ehhfinder.calc_single_ehh(p.query); //query_locus
                }
            }else if(p.CALC_IHS || p.CALC_NSL){
                IHS ihsfinder(hm, p, &flog, &fout);
                ihsfinder.main(); //thread pool
            }else if (p.CALC_XP || p.CALC_XPNSL)
            {
                XPIHH xpihhfinder(hm, p, &flog, &fout);
                xpihhfinder.main();
            }else if (p.CALC_SOFT){
                IHH12 ihh12finder(hm, p, &flog, &fout);
                ihh12finder.main();
            }else if(p.CALC_PI){
                PI pifinder(hm, p, &flog, &fout);
                pifinder.main();
            }
            flog.close();
            fout.close();

        // #ifdef PTW32_STATIC_LIB
        //     pthread_win32_process_detach_np();
        // #endif

            double time_end = (clock() / (double) CLOCKS_PER_SEC);
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = end - start;
            std::cout<<("Total CPU time: "+to_string(time_end-time_start)+" s. \n");
            std::cout << "Program took " << duration.count() << " seconds to complete." << std::endl;

            return 0;
        }
};

#endif