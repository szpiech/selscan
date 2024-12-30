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

            p.LOW_MEM = true;
            //params.getBoolFlag(ARG_LOW_MEM);

            p.MISSING = params.getBoolFlag(ARG_MISSING_FLAG);

            
            p.MAX_EXTEND_NSL = params.getIntFlag(ARG_MAX_EXTEND_NSL); //p.MAX_EXTEND_NSL = ( params.getIntFlag(ARG_MAX_EXTEND_NSL) <= 0 ) ? hm.mapData.nloci : params.getIntFlag(ARG_MAX_EXTEND_NSL);
            p.MAX_EXTEND =  params.getIntFlag(ARG_MAX_EXTEND);
            p.TRUNC = params.getBoolFlag(ARG_TRUNC);
            p.CHR_LIST = params.getStringFlag(ARG_CHR_LIST);

            
            // this->hm = hm;
            //this->bar = new Bar();
            //this->numThreads = p.numThreads
        }
        
        ~MainTools(){
            //delete hm;
        }

        int run(int argc, char *argv[])
        {
            auto start = std::chrono::high_resolution_clock::now();
            double time_start = (clock() / (double)CLOCKS_PER_SEC);

            char PI_WIN_str[50];
            snprintf(PI_WIN_str, 50, "%d", p.PI_WIN);

            if (p.query.compare(DEFAULT_EHH) != 0)
                p.SINGLE_EHH = true;

            // input data is loaded into HapMap object
            // hm = new HapMap(p, &flog, &fout);
            // std::unique_ptr<HapMap>

            hm = std::make_unique<HapMap>(p);
            bool ERROR = hm->loadHapMapData();
            if (ERROR)
            {
                cerr << "ERROR: could not load hapmap data.\n";
                throw 1;
            }

            for (int i = 0; i < argc; i++)
            {
                // flog << argv[i] << " ";
                cerr << argv[i] << " ";
            }

            // flog << "\nv" + VERSION + "\nCalculating ";
            cerr << "\nv" + VERSION + "\nCalculating ";
            if (p.CALC_XP)
                cerr << "XP-EHH.\n";
            if (p.CALC_PI)
                cerr << "PI.\n";
            if (p.CALC_IHS)
                cerr << " iHS.\n";
            if (p.CALC_NSL)
                cerr << " nSL.\n";
            if (p.CALC_XPNSL)
                cerr << " XP-nSL.\n";
            if (p.CALC_SOFT)
                cerr << " iHH1K.\n";

            if (hm->mapData->nloci < p.numThreads)
            {
                p.numThreads = 1;
                cerr << "WARNING: there are fewer loci than threads requested.  Running with " << p.numThreads << " thread instead.\n";
            }

            /*
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
            */

            // @@ EXPERIMENTAL: instead of else if allow all to run
            if (p.SINGLE_EHH)
            {
                EHH ehhfinder(hm, p);
                if (p.CALC_SOFT)
                {
                    throw("ERROR: Soft EHH not implemented yet.\n");
                }
                else
                {
                    ehhfinder.calc_single_ehh(p.query); // query_locus
                }
            }
            if (p.CALC_IHS)
            {
                bool stored_calc_nsl = p.CALC_NSL;
                p.CALC_NSL = false;
                IHS ihsfinder(hm, p);
                ihsfinder.main(); // thread pool
                p.CALC_NSL = stored_calc_nsl;
                p.CALC_IHS = false;
            }
            if (p.CALC_NSL)
            {
                IHS ihsfinder(hm, p);
                ihsfinder.main(); // thread pool
            }

            if (p.CALC_XP)
            {
                bool stored_calc_xpnsl = p.CALC_XPNSL;
                p.CALC_XPNSL = false;
                XPIHH xpihhfinder(hm, p);
                xpihhfinder.main();
                p.CALC_XPNSL = stored_calc_xpnsl;
                p.CALC_XP = false;
            }

            if (p.CALC_XPNSL)
            {
                XPIHH xpihhfinder(hm, p);
                xpihhfinder.main();
            }

            if (p.CALC_SOFT)
            {
                IHH12 ihh12finder(hm, p);
                ihh12finder.main();
            }

            if (p.CALC_PI)
            {
                PI pifinder(hm, p);
                pifinder.main();
            }

            // #ifdef PTW32_STATIC_LIB
            //     pthread_win32_process_detach_np();
            // #endif

            double time_end = (clock() / (double)CLOCKS_PER_SEC);
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = end - start;
            std::cout << ("Total CPU time: " + to_string(time_end - time_start) + " s. \n");
            std::cout << "Program took " << duration.count() << " seconds to complete." << std::endl;

            return 0;
            
        }
};

#endif