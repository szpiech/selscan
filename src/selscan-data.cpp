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
#include "selscan-data.h"
#include "selscan-cli.h"

// #include "selscan-maintools.h"
#include <sstream>
#include <queue> 
#include<algorithm>


using namespace std;


//Returns 1 on error
bool HapMap::loadHapMapData(param_main &p, int argc, char *argv[], ofstream* flog, ofstream* fout){
    this->flog = flog;
    hapData.flog = flog;
    hapData2.flog = flog;


    this->fout = fout;

    string hapFilename = p.hapFilename;
    string hapFilename2 =  p.hapFilename2;
    string mapFilename = p.mapFilename;
    string tpedFilename = p.tpedFilename;
    string tpedFilename2 = p.tpedFilename2;
    string vcfFilename = p.vcfFilename;
    string vcfFilename2 = p.vcfFilename2;
    
    bool TPED = false;
    if (tpedFilename.compare(DEFAULT_FILENAME_POP1_TPED) != 0) TPED = true;

    bool VCF = false;
    if (vcfFilename.compare(DEFAULT_FILENAME_POP1_VCF) != 0) VCF = true;

    bool UNPHASED = p.UNPHASED;
    
    bool USE_PMAP = p.USE_PMAP;
    bool ALT = p.ALT;
    bool WAGH = p.WAGH;
    bool CALC_IHS = p.CALC_IHS;
    bool CALC_XPNSL = p.CALC_XPNSL;
    bool CALC_NSL = p.CALC_NSL;
    bool WRITE_DETAILED_IHS = p.WRITE_DETAILED_IHS;
    bool CALC_XP = p.CALC_XP;
    bool CALC_SOFT = p.CALC_SOFT;
    bool SINGLE_EHH = p.SINGLE_EHH;
    bool LOW_MEM = p.LOW_MEM;

    
    hapData.initParams(UNPHASED, p.SKIP, p.MAF, p.numThreads);
    hapData2.initParams(UNPHASED, p.SKIP, p.MAF, p.numThreads);

    auto start_reading = std::chrono::high_resolution_clock::now();

    if(LOW_MEM){
        try
        {
            if (VCF) {
                hapData.readHapDataVCF_bitset(vcfFilename);
                if (CALC_XP || CALC_XPNSL)
                {
                    hapData2.readHapDataVCF_bitset(vcfFilename2);
                    if (hapData.nloci != hapData2.nloci)
                    {
                        std::cerr << "ERROR: Haplotypes from " << vcfFilename << " and " << vcfFilename2 << " do not have the same number of loci.\n";
                        return 1;
                    }
                }
                if(!CALC_NSL && !CALC_XPNSL && !USE_PMAP) {
                    mapData.readMapData(mapFilename, hapData.nloci, USE_PMAP, hapData.skipQueue);
                }
                else{//Load physical positions
                    mapData.readMapDataVCF(vcfFilename, hapData.nloci, hapData.skipQueue);
                }
            }
            else
            {
                if (CALC_XP || CALC_XPNSL)
                {
                    hapData.initParams(UNPHASED, false, p.MAF, p.numThreads);
                    hapData2.initParams(UNPHASED, false, p.MAF, p.numThreads);

                    hapData.readHapData_bitset(hapFilename);
                    hapData2.readHapData_bitset(hapFilename2);
                    if (hapData.nloci != hapData2.nloci)
                    {
                        std::cerr << "ERROR: Haplotypes from " << hapFilename << " and " << hapFilename2 << " do not have the same number of loci.\n";
                        return 1;
                    }
                }else{
                    hapData.readHapData_bitset(hapFilename);
                }
                mapData.readMapData(mapFilename, hapData.nloci, USE_PMAP, hapData.skipQueue);
                
            }
        }
        catch (...)
        {
            return 1;
        }
    }else{
        try
        {
            if (TPED)
            {
                hapData.readHapDataTPED(tpedFilename);
                if (CALC_XP || CALC_XPNSL)
                {
                    hapData2.readHapDataTPED(tpedFilename2);
                    if (hapData.nloci != hapData2.nloci)
                    {
                        std::cerr << "ERROR: Haplotypes from " << tpedFilename << " and " << tpedFilename2 << " do not have the same number of loci.\n";
                        return 1;
                    }
                }
                mapData.readMapDataTPED(tpedFilename, hapData.nloci, hapData.nhaps, USE_PMAP);
            }
            else if (VCF) {
                handleData(vcfFilename, "VCF");
                
                //exit(1);
                //ataReader dr(vcfFilename, hapData);
                
                
                //hapData.readHapDataVCF(vcfFilename);
                
                if (CALC_XP || CALC_XPNSL)
                {
                    hapData2.readHapDataVCF(vcfFilename2);
                    if (hapData.nloci != hapData2.nloci)
                    {
                        std::cerr << "ERROR: Haplotypes from " << vcfFilename << " and " << vcfFilename2 << " do not have the same number of loci.\n";
                        return 1;
                    }
                }
                if(!CALC_NSL && !CALC_XPNSL && !USE_PMAP) {
                    mapData.readMapData(mapFilename, hapData.nloci, USE_PMAP, hapData.skipQueue);
                }
                else{//Load physical positions
                    mapData.readMapDataVCF(vcfFilename, hapData.nloci, hapData.skipQueue);
                }
            }
            else
            {
                if (CALC_XP || CALC_XPNSL)
                {
                    hapData.initParams(UNPHASED, false, p.MAF, p.numThreads);
                    hapData2.initParams(UNPHASED, false, p.MAF, p.numThreads);

                    hapData.readHapData(hapFilename);
                    hapData2.readHapData(hapFilename2);
                    if (hapData.nloci != hapData2.nloci)
                    {
                        std::cerr << "ERROR: Haplotypes from " << hapFilename << " and " << hapFilename2 << " do not have the same number of loci.\n";
                        return 1;
                    }
                }else{
                    //hapData.readHapData(hapFilename);
                    handleData(hapFilename, "HAP");
                }
                mapData.readMapData(mapFilename, hapData.nloci, USE_PMAP, hapData.skipQueue);
                
            }
        }
        catch (...)
        {
            return 1;
        }
    }

    auto end_reading = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> read_duration =  end_reading - start_reading;
    //FLOG
    cout<<("Input file loaded in "+to_string(read_duration.count())+" s.\n");

    (*(flog))<<("Input file loaded in "+to_string(read_duration.count())+" s.\n");
    //mapData.print();


    // Check if map is in order
    for (int i = 1; i < mapData.nloci; i++) {
        if ( mapData.mapEntries[i].physicalPos <  mapData.mapEntries[i-1].physicalPos ) {
            std::cerr << "ERROR: Variant physical position must be monotonically increasing.\n";
            std::cerr << "\t" << mapData.mapEntries[i].locusName << " " << mapData.mapEntries[i].physicalPos << " appears after";
            std::cerr << "\t" <<  mapData.mapEntries[i-1].locusName << " " << mapData.mapEntries[i-1].physicalPos << "\n";
            return 1;
        }
        if ( !CALC_NSL && mapData.mapEntries[i].geneticPos  < mapData.mapEntries[i-1].geneticPos  ) {
            std::cerr << "ERROR: Variant genetic position must be monotonically increasing.\n";
            std::cerr << "\t" << mapData.mapEntries[i].locusName << " " << mapData.mapEntries[i].geneticPos << " appears after";
            std::cerr << "\t" << mapData.mapEntries[i-1].locusName << " " << mapData.mapEntries[i-1].geneticPos << "\n";
            return 1;
        }
    }

    return 0;
}
