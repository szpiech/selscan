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
#include "selscan-maintools.h"
#include "selscan-data.h"
#include "selscan-pbar.h"
#include "param_t.h"

using namespace std;

int main(int argc, char *argv[])
{
    cerr << "selscan v" + VERSION + "\n";
#ifdef PTW32_STATIC_LIB
    pthread_win32_process_attach_np();
#endif

    param_t params;
    
    bool ERROR = initalizeParameters(params,argc,argv);

    if (ERROR) return 1;

    ERROR = checkParameters(params,argc,argv);

    if (ERROR) return 1;

    string hapFilename = params.getStringFlag(ARG_FILENAME_POP1);
    string hapFilename2 = params.getStringFlag(ARG_FILENAME_POP2);
    string mapFilename = params.getStringFlag(ARG_FILENAME_MAP);
    string tpedFilename = params.getStringFlag(ARG_FILENAME_POP1_TPED);
    string tpedFilename2 = params.getStringFlag(ARG_FILENAME_POP2_TPED);
    string vcfFilename = params.getStringFlag(ARG_FILENAME_POP1_VCF);
    string vcfFilename2 = params.getStringFlag(ARG_FILENAME_POP2_VCF);
    
    bool TPED = false;
    if (tpedFilename.compare(DEFAULT_FILENAME_POP1_TPED) != 0) TPED = true;

    bool VCF = false;
    if (vcfFilename.compare(DEFAULT_FILENAME_POP1_VCF) != 0) VCF = true;

    string outFilename = params.getStringFlag(ARG_OUTFILE);
    string query = params.getStringFlag(ARG_EHH);

    int queryLoc = -1;
    int numThreads = params.getIntFlag(ARG_THREAD);
    int SCALE_PARAMETER = params.getIntFlag(ARG_GAP_SCALE);
    int MAX_GAP = params.getIntFlag(ARG_MAX_GAP);

    double EHH_CUTOFF = params.getDoubleFlag(ARG_CUTOFF);
    double MAF = params.getDoubleFlag(ARG_MAF);

    bool UNPHASED = params.getBoolFlag(ARG_UNPHASED);
    bool USE_PMAP = params.getBoolFlag(ARG_PMAP);
    bool ALT = params.getBoolFlag(ARG_ALT);
    bool WAGH = params.getBoolFlag(ARG_WAGH);
    bool CALC_IHS = params.getBoolFlag(ARG_IHS);
    bool CALC_XPNSL = params.getBoolFlag(ARG_XPNSL);
    bool CALC_NSL = params.getBoolFlag(ARG_NSL);
    bool WRITE_DETAILED_IHS = params.getBoolFlag(ARG_IHS_DETAILED);
    bool CALC_XP = params.getBoolFlag(ARG_XP);
    bool CALC_SOFT = params.getBoolFlag(ARG_SOFT);
    bool SINGLE_EHH = false;
    
    bool SKIP = !params.getBoolFlag(ARG_KEEP);//params.getBoolFlag(ARG_SKIP);
    if(params.getBoolFlag(ARG_SKIP)){
        cerr << "WARNING: " << ARG_SKIP << " is now on by dafault.  This flag no longer has a function.\n";
    }
    //bool TRUNC = params.getBoolFlag(ARG_TRUNC);

    int EHH1K = params.getIntFlag(ARG_SOFT_K);

    bool CALC_PI = params.getBoolFlag(ARG_PI);
    int PI_WIN = params.getIntFlag(ARG_PI_WIN);

    char PI_WIN_str[50];
    snprintf(PI_WIN_str,50, "%d", PI_WIN);

    if (query.compare(DEFAULT_EHH) != 0) SINGLE_EHH = true;


    if (SINGLE_EHH) outFilename += ".ehh." + query;
    else if (CALC_IHS) outFilename += ".ihs";
    else if (CALC_NSL) outFilename += ".nsl";
    else if (CALC_XPNSL) outFilename += ".xpnsl";
    else if (CALC_XP) outFilename += ".xpehh";
    else if (CALC_SOFT) outFilename += ".ihh12";
    else if (CALC_PI) outFilename += ".pi." + string(PI_WIN_str) + "bp";

    if (ALT) outFilename += ".alt";

    HaplotypeData *hapData, *hapData2;
    MapData *mapData;

    ERROR = loadHapMapData(&hapData,&hapData2,&mapData,params,argc,argv);

    if (ERROR) return 1;

    //Open stream for log file
    ofstream flog;
    string logFilename = outFilename + ".log";
    flog.open(logFilename.c_str());
    if (flog.fail())
    {
        cerr << "ERROR: could not open " << logFilename << " for writing.\n";
        return 1;
    }

    //Open stream for output file
    ofstream fout;
    outFilename += ".out";
    fout.open(outFilename.c_str());
    if (fout.fail())
    {
        cerr << "ERROR: could not open " << outFilename << " for writing.\n";
        return 1;
    }


    for (int i = 0; i < argc; i++)
    {
        flog << argv[i] << " ";
    }
    flog << "\nv" + VERSION + "\nCalculating ";
    if (CALC_XP) flog << "XP-EHH.\n";
    else if (CALC_PI) flog << "PI.\n";
    else if (CALC_IHS) flog << " iHS.\n";
    else if (CALC_NSL) flog << " nSL.\n";
    else if (CALC_XPNSL) flog << " XP-nSL.\n";
    else if (CALC_SOFT) flog << " iHH1K.\n";

    if(params.getBoolFlag(ARG_SKIP)){
        flog << "WARNING: " << ARG_SKIP << " is now on by dafault.  This flag no longer has a function.\n";
    }

    if (TPED)
    {
        flog << "Input filename: " << tpedFilename << "\n";
        if (CALC_XP || CALC_XPNSL) flog << "Reference input filename: " << tpedFilename2 << "\n";

    }
    else if (VCF) {
        flog << "Input filename: " << vcfFilename << "\n";
        if (CALC_XP || CALC_XPNSL) flog << "Reference input filename: " << vcfFilename2 << "\n";
        flog << "Map filename: " << mapFilename << "\n";
    }
    else {
        flog << "Input filename: " << hapFilename << "\n";
        if (CALC_XP || CALC_XPNSL) flog << "Reference input filename: " << hapFilename2 << "\n";
        flog << "Map filename: " << mapFilename << "\n";
    }
    flog << "Output file: " << outFilename << "\n";
    flog << "Threads: " << numThreads << "\n";
    flog << "Scale parameter: " << SCALE_PARAMETER << "\n";
    flog << "Max gap parameter: " << MAX_GAP << "\n";
    flog << "EHH cutoff value: " << EHH_CUTOFF << "\n";
    flog << "Phased: ";
    if(UNPHASED) flog << "no\n";
    else flog << "yes\n";
    flog << "Alt flag set: ";
    if (ALT) flog << "yes\n";
    else flog << "no\n";
    flog.flush();

    Bar pbar;

    double *ihs, *ihh1, *ihh2;
    double *ihhDerivedLeft, *ihhDerivedRight, *ihhAncestralLeft, *ihhAncestralRight;
    double *freq, *freq1, *freq2;

    if (mapData->nloci < numThreads)
    {
        numThreads = 1;
        cerr << "WARNING: there are fewer loci than threads requested.  Running with " << numThreads << " thread instead.\n";
    }

    if (SINGLE_EHH)
    {

        freq = new double[hapData->nloci];

        MapData *newMapData;
        HaplotypeData *newHapData;
        double *newfreq;

        int count = 0;
        for (int i = 0; i < hapData->nloci; i++)
        {
            freq[i] = calcFreq(hapData, i, UNPHASED);
            if (freq[i] > MAF && 1 - freq[i] > MAF) count++;
        }


        queryLoc = queryFound(mapData, query);
        double queryFreq = calcFreq(hapData, queryLoc, UNPHASED);
        if (queryLoc < 0)
        {
            cerr << "ERROR: Could not find specific locus query, " << query << ", in data.\n";
            return 1;
        }
        else if (SKIP && (queryFreq < MAF || 1 - queryFreq < MAF))
        {
            cerr << "ERROR: EHH not calculated for " << query << ". MAF < " << MAF << ".\n";
            cerr << "\tIf you wish to calculate EHH at this locus either change --maf or set --keep-low-freq.\n";
            return 1;
        }
        else if (!SKIP && (queryFreq == 0 || queryFreq == 1)){
            cerr << "ERROR: EHH not calculated for " << query << ". Frequency = " << queryFreq << ".\n";
            return 1;
        }
        else
        {
            cerr << "Found " << query << " in data.\n";
        }

        if (SKIP) //prefilter all sites < MAF
        {
            cerr << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
            flog << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
            newfreq = new double [count];
            newMapData = initMapData(count);
            newMapData->chr = mapData->chr;
            int j = 0;
            for (int locus = 0; locus < mapData->nloci; locus++)
            {
                if (freq[locus] <= MAF || 1 - freq[locus] <= MAF)
                {
                    continue;
                }
                else
                {
                    newMapData->physicalPos[j] = mapData->physicalPos[locus];
                    newMapData->geneticPos[j] = mapData->geneticPos[locus];
                    newMapData->locusName[j] = mapData->locusName[locus];
                    newfreq[j] = freq[locus];
                    j++;
                }
            }

            newHapData = initHaplotypeData(hapData->nhaps, count);

            for (int hap = 0; hap < newHapData->nhaps; hap++)
            {
                j = 0;
                for (int locus = 0; locus < mapData->nloci; locus++)
                {
                    if (freq[locus] <= MAF || 1 - freq[locus] <= MAF)
                    {
                        continue;
                    }
                    else
                    {
                        newHapData->data[hap][j] = hapData->data[hap][locus];
                        j++;
                    }
                }
            }

            cerr << "Removed " << mapData->nloci - count << " low frequency variants.\n";
            flog << "Removed " << mapData->nloci - count << " low frequency variants.\n";

            delete [] freq;
            freq = newfreq;
            newfreq = NULL;

            releaseHapData(hapData);
            hapData = newHapData;
            newHapData = NULL;

            releaseMapData(mapData);
            mapData = newMapData;
            newMapData = NULL;
        }

        queryLoc = queryFound(mapData, query);

        work_order_t *order = new work_order_t;
        pthread_t *peer = new pthread_t;
        order->hapData = hapData;
        order->mapData = mapData;
        order->flog = &flog;
        order->fout = &fout;
        order->filename = outFilename;
        order->params = &params;
        order->queryLoc = queryLoc;
        order->calc = &calculateHomozygosity;

        if (CALC_SOFT)
        {
            pthread_create(peer,
                           NULL,
                           (void *(*)(void *))query_locus_soft,
                           (void *)order);
            pthread_join(*peer, NULL);
        }
        else
        {
            pthread_create(peer,
                           NULL,
                           (void *(*)(void *))query_locus,
                           (void *)order);
            pthread_join(*peer, NULL);
        }

        delete peer;
        delete order;
        return 0;
    }

    if (CALC_XP || CALC_XPNSL)
    {

        freq1 = new double[hapData->nloci];
        freq2 = new double[hapData2->nloci];

        for (int i = 0; i < hapData->nloci; i++)
        {
            freq1[i] = calcFreq(hapData, i, UNPHASED);
            freq2[i] = calcFreq(hapData2, i, UNPHASED);
        }


        ihh1 = new double[mapData->nloci];
        ihh2 = new double[mapData->nloci];
        
        barInit(pbar, mapData->nloci, 78);

        if (CALC_XPNSL){
            for (int i = 0; i < mapData->nloci; i++){
                mapData->geneticPos[i] = i;
            }
        }

        if (CALC_XP) cerr << "Starting XP-EHH calculations.\n";
        if (CALC_XPNSL) cerr << "Starting XP-nSL calculations.\n";
        work_order_t *order;
        pthread_t *peer = new pthread_t[numThreads];
        //int prev_index = 0;
        for (int i = 0; i < numThreads; i++)
        {
            order = new work_order_t;
            order->id = i;
            order->hapData1 = hapData;
            order->hapData2 = hapData2;
            order->mapData = mapData;
            order->ihh1 = ihh1;
            order->ihh2 = ihh2;
            order->freq1 = freq1;
            order->freq2 = freq2;
            order->flog = &flog;
            order->bar = &pbar;
            order->params = &params;
            pthread_create(&(peer[i]),
                           NULL,
                           (void *(*)(void *))calc_xpihh,
                           (void *)order);
        }

        for (int i = 0; i < numThreads; i++)
        {
            pthread_join(peer[i], NULL);
        }

        delete [] peer;
        releaseHapData(hapData);
        releaseHapData(hapData2);
        cerr << "\nFinished.\n";

        if (CALC_XP) fout << "id\tpos\tgpos\tp1\tihh1\tp2\tihh2\txpehh\n";
        if (CALC_XPNSL) fout << "id\tpos\tgpos\tp1\tsL1\tp2\tsL2\txpnsl\n";
        for (int i = 0; i < mapData->nloci; i++)
        {
            if (ihh1[i] != MISSING && ihh2[i] != MISSING && ihh1[i] != 0 && ihh2[i] != 0)
            {
                fout << mapData->locusName[i] << "\t"
                     << mapData->physicalPos[i] << "\t"
                     << mapData->geneticPos[i] << "\t"
                     << freq1[i] << "\t"
                     << ihh1[i] << "\t"
                     << freq2[i] << "\t"
                     << ihh2[i] << "\t";
                fout << log10(ihh1[i] / ihh2[i]) << endl;
            }
        }
    }
    else if (CALC_IHS)
    {

        freq = new double[hapData->nloci];

        MapData *newMapData;
        HaplotypeData *newHapData;
        double *newfreq;

        int count = 0;
        for (int i = 0; i < hapData->nloci; i++)
        {
            freq[i] = calcFreq(hapData, i, UNPHASED);
            if (freq[i] > MAF && 1 - freq[i] > MAF) count++;
        }

        if (SKIP) //prefilter all sites < MAF
        {
            cerr << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
            flog << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
            newfreq = new double [count];
            newMapData = initMapData(count);
            newMapData->chr = mapData->chr;
            int j = 0;
            for (int locus = 0; locus < mapData->nloci; locus++)
            {
                if (freq[locus] <= MAF || 1 - freq[locus] <= MAF)
                {
                    continue;
                }
                else
                {
                    newMapData->physicalPos[j] = mapData->physicalPos[locus];
                    newMapData->geneticPos[j] = mapData->geneticPos[locus];
                    newMapData->locusName[j] = mapData->locusName[locus];
                    newfreq[j] = freq[locus];
                    j++;
                }
            }

            newHapData = initHaplotypeData(hapData->nhaps, count);

            for (int hap = 0; hap < newHapData->nhaps; hap++)
            {
                j = 0;
                for (int locus = 0; locus < mapData->nloci; locus++)
                {
                    if (freq[locus] <= MAF || 1 - freq[locus] <= MAF)
                    {
                        continue;
                    }
                    else
                    {
                        newHapData->data[hap][j] = hapData->data[hap][locus];
                        j++;
                    }
                }
            }

            cerr << "Removed " << mapData->nloci - count << " low frequency variants.\n";
            flog << "Removed " << mapData->nloci - count << " low frequency variants.\n";

            delete [] freq;
            freq = newfreq;
            newfreq = NULL;

            releaseHapData(hapData);
            hapData = newHapData;
            newHapData = NULL;

            releaseMapData(mapData);
            mapData = newMapData;
            newMapData = NULL;
        }

        ihh1 = new double[mapData->nloci];
        ihh2 = new double[mapData->nloci];
        ihs = new double[hapData->nloci];

        if (WRITE_DETAILED_IHS) {
            ihhDerivedLeft = new double[hapData->nloci];
            ihhDerivedRight = new double[hapData->nloci];
            ihhAncestralLeft = new double[hapData->nloci];
            ihhAncestralRight = new double[hapData->nloci];
        }
        barInit(pbar, mapData->nloci, 78);

        cerr << "Starting iHS calculations with alt flag ";
        if (!ALT) cerr << "not ";
        cerr << "set.\n";

        work_order_t *order;
        pthread_t *peer = new pthread_t[numThreads];
        //int prev_index = 0;
        for (int i = 0; i < numThreads; i++)
        {
            order = new work_order_t;
            order->id                = i;
            order->hapData           = hapData;
            order->mapData           = mapData;
            order->ihh1              = ihh1;
            order->ihh2              = ihh2;
            order->ihs               = ihs;
            order->ihhDerivedLeft    = ihhDerivedLeft;
            order->ihhDerivedRight   = ihhDerivedRight;
            order->ihhAncestralLeft  = ihhAncestralLeft;
            order->ihhAncestralRight = ihhAncestralRight;
            order->freq              = freq;
            order->flog              = &flog;
            order->bar               = &pbar;
            order->params            = &params;
            order->calc              = &calculateHomozygosity;

            pthread_create(&(peer[i]),
                           NULL,
                           (void *(*)(void *))calc_ihs,
                           (void *)order);
        }

        for (int i = 0; i < numThreads; i++)
        {
            pthread_join(peer[i], NULL);
        }

        delete [] peer;
        releaseHapData(hapData);
        cerr << "\nFinished.\n";

        for (int i = 0; i < mapData->nloci; i++)
        {
            if (ihs[i] != MISSING && ihh1[i] != 0 && ihh2[i] != 0)
            {
                fout << mapData->locusName[i] << "\t"
                     << mapData->physicalPos[i] << "\t"
                     << freq[i] << "\t"
                     << ihh1[i] << "\t"
                     << ihh2[i] << "\t"
                     << ihs[i];
                if (!WRITE_DETAILED_IHS)
                {
                    fout << endl;
                } else
                {
                    fout << "\t"
                         << ihhDerivedLeft[i]    << "\t"
                         << ihhDerivedRight[i]   << "\t"
                         << ihhAncestralLeft[i]  << "\t"
                         << ihhAncestralRight[i] << endl;
                }
            }
        }
    }
    else if (CALC_NSL)
    {

        freq = new double[hapData->nloci];

        MapData *newMapData;
        HaplotypeData *newHapData;
        double *newfreq;

        int count = 0;
        for (int i = 0; i < hapData->nloci; i++)
        {
            freq[i] = calcFreq(hapData, i, UNPHASED);
            if (freq[i] > MAF && 1 - freq[i] > MAF) count++;
        }
        
        if (SKIP) //prefilter all sites < MAF
        {
            cerr << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
            flog << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
            newfreq = new double [count];
            newMapData = initMapData(count);
            newMapData->chr = mapData->chr;
            int j = 0;
            for (int locus = 0; locus < mapData->nloci; locus++)
            {
                if (freq[locus] <= MAF || 1 - freq[locus] <= MAF)
                {
                    continue;
                }
                else
                {
                    newMapData->physicalPos[j] = mapData->physicalPos[locus];
                    newMapData->geneticPos[j] = mapData->geneticPos[locus];
                    newMapData->locusName[j] = mapData->locusName[locus];
                    newfreq[j] = freq[locus];
                    j++;
                }
            }
            newHapData = initHaplotypeData(hapData->nhaps, count);
            for (int hap = 0; hap < newHapData->nhaps; hap++)
            {
                j = 0;
                for (int locus = 0; locus < mapData->nloci; locus++)
                {
                    if (freq[locus] <= MAF || 1 - freq[locus] <= MAF)
                    {
                        continue;
                    }
                    else
                    {
                        newHapData->data[hap][j] = hapData->data[hap][locus];
                        j++;
                    }
                }
            }
            cerr << "Removed " << mapData->nloci - count << " low frequency variants.\n";
            flog << "Removed " << mapData->nloci - count << " low frequency variants.\n";
            delete [] freq;
            freq = newfreq;
            newfreq = NULL;
            releaseHapData(hapData);
            hapData = newHapData;
            newHapData = NULL;
            releaseMapData(mapData);
            mapData = newMapData;
            newMapData = NULL;
        }

        ihh1 = new double[mapData->nloci];
        ihh2 = new double[mapData->nloci];
        ihs = new double[hapData->nloci];

        barInit(pbar, mapData->nloci, 78);

        cerr << "Starting nSL calculations with alt flag ";
        if (!ALT) cerr << "not ";
        cerr << "set.\n";

        work_order_t *order;
        pthread_t *peer = new pthread_t[numThreads];
        //int prev_index = 0;
        for (int i = 0; i < numThreads; i++)
        {
            order = new work_order_t;
            order->id = i;
            order->hapData = hapData;
            order->mapData = mapData;
            order->ihh1 = ihh1;
            order->ihh2 = ihh2;
            order->ihs = ihs;
            order->freq = freq;
            order->flog = &flog;
            order->bar = &pbar;
            order->params = &params;
            order->calc = &calculateHomozygosity;

            pthread_create(&(peer[i]),
                           NULL,
                           (void *(*)(void *))calc_nsl,
                           (void *)order);
        }

        for (int i = 0; i < numThreads; i++)
        {
            pthread_join(peer[i], NULL);
        }

        delete [] peer;
        releaseHapData(hapData);
        cerr << "\nFinished.\n";

        for (int i = 0; i < mapData->nloci; i++)
        {
            if (ihs[i] != MISSING && ihh1[i] != 0 && ihh2[i] != 0)
            {
                fout << mapData->locusName[i] << "\t"
                     << mapData->physicalPos[i] << "\t"
                     << freq[i] << "\t"
                     << ihh1[i] << "\t"
                     << ihh2[i] << "\t"
                     << ihs[i] << endl;
            }
        }
    }
    else if (CALC_SOFT)
    {
        freq = new double[hapData->nloci];

        MapData *newMapData;
        HaplotypeData *newHapData;
        double *newfreq;

        int count = 0;
        for (int i = 0; i < hapData->nloci; i++)
        {
            freq[i] = calcFreq(hapData, i, UNPHASED);
            if (freq[i] > MAF && 1 - freq[i] > MAF) count++;
        }

        if (SKIP) //prefilter all sites < MAF
        {
            cerr << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
            flog << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
            newfreq = new double [count];
            newMapData = initMapData(count);
            newMapData->chr = mapData->chr;
            int j = 0;
            for (int locus = 0; locus < mapData->nloci; locus++)
            {
                if (freq[locus] <= MAF || 1 - freq[locus] <= MAF)
                {
                    continue;
                }
                else
                {
                    newMapData->physicalPos[j] = mapData->physicalPos[locus];
                    newMapData->geneticPos[j] = mapData->geneticPos[locus];
                    newMapData->locusName[j] = mapData->locusName[locus];
                    newfreq[j] = freq[locus];
                    j++;
                }
            }

            newHapData = initHaplotypeData(hapData->nhaps, count);

            for (int hap = 0; hap < newHapData->nhaps; hap++)
            {
                j = 0;
                for (int locus = 0; locus < mapData->nloci; locus++)
                {
                    if (freq[locus] <= MAF || 1 - freq[locus] <= MAF)
                    {
                        continue;
                    }
                    else
                    {
                        newHapData->data[hap][j] = hapData->data[hap][locus];
                        j++;
                    }
                }
            }

            cerr << "Removed " << mapData->nloci - count << " low frequency variants.\n";
            flog << "Removed " << mapData->nloci - count << " low frequency variants.\n";

            delete [] freq;
            freq = newfreq;
            newfreq = NULL;

            releaseHapData(hapData);
            hapData = newHapData;
            newHapData = NULL;

            releaseMapData(mapData);
            mapData = newMapData;
            newMapData = NULL;
        }

        ihh1 = new double[mapData->nloci];
        ihh2 = new double[mapData->nloci];
        ihs = new double[hapData->nloci];

        barInit(pbar, mapData->nloci, 78);

        cerr << "Starting iHH12 calculations with alt flag ";
        if (!ALT) cerr << "not ";
        cerr << "set.\n";

        work_order_t *order;
        pthread_t *peer = new pthread_t[numThreads];
        //int prev_index = 0;
        for (int i = 0; i < numThreads; i++)
        {
            order = new work_order_t;
            //order->first_index = prev_index;
            //order->last_index = prev_index + NUM_PER_THREAD[i];
            //prev_index += NUM_PER_THREAD[i];
            order->hapData = hapData;
            order->mapData = mapData;
            order->ihh1 = ihh1;
            order->ihh2 = ihh2;
            order->ihs = ihs;
            order->freq = freq;
            order->flog = &flog;
            order->bar = &pbar;
            order->params = &params;
            order->id = i;

            pthread_create(&(peer[i]),
                           NULL,
                           (void *(*)(void *))calc_soft_ihs,
                           (void *)order);
        }

        for (int i = 0; i < numThreads; i++)
        {
            pthread_join(peer[i], NULL);
        }

        delete [] peer;
        releaseHapData(hapData);
        cerr << "\nFinished.\n";

        fout << "id\tpos\tp1\tihh12\n";
        
        for (int i = 0; i < mapData->nloci; i++)
        {
            if (ihs[i] != MISSING )//&& ihh1[i] != MISSING && ihh2[i] != MISSING)
            {
                fout << mapData->locusName[i] << "\t";
                fout << mapData->physicalPos[i] << "\t";
                fout << freq[i] << "\t";
                //fout << ihh1[i] << "\t"; //ihh1
                fout << ihs[i] << "\n"; //ihh12
                //fout << ihh2[i] << endl; //ihh2d1
            }
        }
    }
    else if (CALC_PI)
    {
        //cerr << "Not implemented.\n";
        //return 1;

        cerr << "Starting pi calculations with " << PI_WIN << "bp windows.\n";

        calculatePi(hapData, mapData, PI_WIN, outFilename);

        releaseHapData(hapData);
        cerr << "\nFinished.\n";

    }

    flog.close();
    fout.close();

#ifdef PTW32_STATIC_LIB
    pthread_win32_process_detach_np();
#endif

    return 0;
}
