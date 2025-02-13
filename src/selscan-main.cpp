//selscan-main.cpp
/* selscan -- a program to calculate EHH-based scans for positive selection in genomes
   Copyright (C) 2014  Zachary A Szpiech
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
#include <iostream>
#include <fstream>
#include <string>
#include <cctype>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <pthread.h>
#include <map>
#include "selscan-cli.h"
#include "selscan-data.h"
#include "param_t.h"


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
#include "stats/ehh12.h"


#include <unordered_map>
#include <thread>

#include <ctime>
#include <cmath>

using namespace std;


string getLogFileName(param_main &p){
    string statname = "";
    if(p.SINGLE_EHH){
        statname = "ehh";
    }
    if (p.CALC_XP){
        statname = "xpehh";
    }
    if (p.CALC_PI){
        statname = "pi";
    }
    if (p.CALC_IHS){
        statname = "ihs";
    }
    if (p.CALC_NSL){
        statname = "nsl";
    }
    if (p.CALC_XPNSL){
        statname = "xpnsl";
    }
    if (p.CALC_SOFT){
        statname = "ihh12";
    }
    return p.outFilename + "." + statname + ".log";
}



/// @brief Run the stats requested by the user, here p is the stat specific parameters (different from base p)
void runStat(std::unique_ptr<HapMap>& hm, param_main &p)
{            
    if (hm->mapData->nloci < p.numThreads)
    {
        p.numThreads = hm->mapData->nloci; // p.numThreads = 1;
        cerr << "WARNING: there are fewer loci than threads requested.  Running with " << p.numThreads << " thread instead.\n";
        *(hm->flog) << "WARNING: there are fewer loci than threads requested.  Running with " << p.numThreads << " thread instead.\n";
    }

    int max_supported_thread = std::thread::hardware_concurrency();
    if(p.numThreads > max_supported_thread){
        cerr << "WARNING: Requested "<< p.numThreads<<" threads, but maximum number of threads supported by hardware is "<< max_supported_thread <<"."<<endl;
        *(hm->flog) << "WARNING: Requested "<< p.numThreads<<" threads, but maximum number of threads supported by hardware is "<< max_supported_thread <<"."<<endl;
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
        ehhfinder.main(p.query); // query_locus
    }
    if (p.SINGLE_EHH12){
        EHH12 ehh12finder(hm, p);
        ehh12finder.main(p.query_ehh12); // query_locus
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
}

int main(int argc, char *argv[])
{
    auto start = std::chrono::high_resolution_clock::now();
    double time_start = (clock() / (double)CLOCKS_PER_SEC);
 
    cerr << "selscan v" + VERSION + "\n";
   
    for (int i = 0; i < argc; i++)
    {
        cerr << argv[i] << " ";
    }
    cerr << "\n";

    


// #ifdef PTW32_STATIC_LIB
//     pthread_win32_process_attach_np();
// #endif

    std::unique_ptr<HapMap> hm; 
    ofstream* flog;

    param_t params;
    param_main p;
    try{
        initalizeParameters(params,argc,argv);
        getBaseParamFromCmdLine(params, p);
        checkParameters(p);
    }catch(const std::exception& e){
        cerr<<"ERROR: "<<e.what()<<endl;
        return 1;
    }
    
    ofstream flog_fs;
    flog = &flog_fs;
    // start logging
    string logFilename = getLogFileName(p);
    (*flog).open(logFilename);
    if ((*flog).fail())
    {
        cerr << "ERROR: could not open " << logFilename << " for writing.\n";
        exit(2);
    }

    (*flog) << "selscan v" + VERSION + "\n";
    for (int i = 0; i < argc; i++)
    {
        (*flog) << argv[i] << " ";
    }
    (*flog) << "\n";


    // if(p.MULTI_CHR){
    //     (*flog)<<"WARNING: Running in multi-chromosome mode.\n";
    //     cerr << "WARNING: Running in multi-chromosome mode.\n";
    // }

    if(p.MULTI_PARAMS){
        //TODON
        p.MULTI_MAF = true;
        (*flog)<<"WARNING: Running in multi-parameter mode.\n";
        cerr<<"WARNING: Running in multi-parameter mode.\n";   
        vector<param_main> multi_params;
        jsonParse(p, multi_params);
        hm = std::make_unique<HapMap>(p, multi_params, flog);
        hm->loadHapMapData();

        // run for all param combinations
        cerr<<"DEBUG::: Running for all param combinations. "<<multi_params.size()<<"\n";
        for(int i = 0; i< multi_params.size(); i++){
            param_main p = multi_params[i];
            p.outFilename = p.outFilename + ".pconfig" + to_string(i);
            runStat(hm, p); // run the main program
            (*flog)<<"DEBUG::: Finished "<<p.outFilename<<".\n";
            cerr<<"DEBUG::: Finished "<<p.outFilename<<".\n";
        }
    }else{
        hm = std::make_unique<HapMap>(p, flog);
        hm->loadHapMapData();
        runStat(hm, p); // run the main program
    }

    double time_end = (clock() / (double)CLOCKS_PER_SEC);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cerr << ("Total CPU time: " + to_string(time_end - time_start) + " s. \n");
    std::cerr << "Program took " << duration.count() << " seconds to complete." << std::endl;
    *flog << ("Total CPU time: " + to_string(time_end - time_start) + " s. \n");
    *flog << "Program took " << duration.count() << " seconds to complete." << std::endl;
    
    // // END LOGGING
    flog->close();


    // #ifdef PTW32_STATIC_LIB
    //     pthread_win32_process_detach_np();
    // #endif
    return 0;
}
