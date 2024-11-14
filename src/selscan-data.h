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

#ifndef __DATA_H__
#define __DATA_H__


#include "selscan-cli.h"

#include <sstream>
#include <algorithm>

#include <string>
#include <iostream>
#include <fstream>
#include <omp.h>
#include "gzstream.h"
#include "param_t.h"

#include <chrono>
#include <queue>
#include <cmath>

#include "hapmap/bitset.h"
#include "hapmap/hapdata.h"
#include "hapmap/mapdata.h"

#include <memory>
#include <atomic>

// #include "filetype/vcf.h"
// #include "filetype/vcf_serial.h"
// #include "filetype/hap.h"
// #include "filetype/hap_serial.h"

using namespace std;

const int MISSING = -9999;  // TO_BE_DELETED
const char MISSING_CHAR = '9';  // TO_BE_DELETED

class HapMap{
public:
    param_main p;
    std::unique_ptr<MapData> mapData;
    std::unique_ptr<HapData> hapData;
    std::unique_ptr<HapData> hapData2;

    ofstream* flog;
    ofstream* fout;

    std::atomic<int> currentProcessed = 0;

    

    HapMap(param_main &p, ofstream* flog, ofstream* fout){
        this->p = p;
        this->flog = flog;
        this->fout = fout;
    }

    ~HapMap(){
    }

    /// Returns 1 on error, Returns 0 on success
    bool loadHapMapData(){
        const string& hapFilename = p.hapFilename;
        const string& hapFilename2 =  p.hapFilename2;
        const string& mapFilename = p.mapFilename;
        const string& tpedFilename = p.tpedFilename;
        const string& tpedFilename2 = p.tpedFilename2;
        const string& vcfFilename = p.vcfFilename;
        const string& vcfFilename2 = p.vcfFilename2;
        
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
        // MISSING
        // if(p.MISSING){
        //     hapData->benchmark_flag = "MISSING";
        // }

        
        mapData = std::make_unique<MapData>(); 
        hapData = std::make_unique<HapData>();
        
        if (CALC_XP || CALC_XPNSL){
            hapData2 = std::make_unique<HapData>();
            hapData->initParams(UNPHASED, false, 0, p.numThreads, flog, p.LOW_MEM, p.MISSING);
            hapData2->initParams(UNPHASED, false, 0, p.numThreads, flog, p.LOW_MEM, p.MISSING);
        }   

        auto start_reading = std::chrono::high_resolution_clock::now();

        
        try
        {
            if (TPED)
            {
                hapData->readHapDataTPED(tpedFilename);
                //handleData(tpedFilename, "TPED", *hapData);
                if (CALC_XP || CALC_XPNSL)
                {
                    hapData2->readHapDataTPED(tpedFilename);
                    //handleData(tpedFilename2, "TPED", *hapData2);
                    if (hapData->nloci != hapData2->nloci)
                    {
                        std::cerr << "ERROR: Haplotypes from " << tpedFilename << " and " << tpedFilename2 << " do not have the same number of loci.\n";
                        return 1;
                    }
                }
                mapData->readMapDataTPED(tpedFilename, hapData->nloci, hapData->nhaps, USE_PMAP, hapData->skipQueue);
            }
            else if (VCF) {
                if (CALC_XP || CALC_XPNSL)
                {
                    hapData->initParams(UNPHASED, false, 0, p.numThreads, flog, p.LOW_MEM);
                    hapData2->initParams(UNPHASED, false, 0, p.numThreads, flog, p.LOW_MEM);

                    hapData->readHapDataVCF(vcfFilename);
                    hapData2->readHapDataVCF(vcfFilename2);

                    // handleData(vcfFilename, "VCF", *hapData);
                    // handleData(vcfFilename2, "VCF", *hapData2);


                    if (hapData->nloci != hapData2->nloci)
                    {
                        std::cerr << "ERROR: Haplotypes from " << vcfFilename << " and " << vcfFilename2 << " do not have the same number of loci.\n";
                        return 1;
                    }
                }else{
                    hapData->initParams(UNPHASED, p.SKIP, p.MAF, p.numThreads, flog, p.LOW_MEM);
                    hapData->readHapDataVCF(vcfFilename);

                    // handleData(vcfFilename, "VCF", *hapData);

                }
                if(!CALC_NSL && !CALC_XPNSL && !USE_PMAP) {
                    mapData->readMapData(mapFilename, hapData->nloci, USE_PMAP, hapData->skipQueue);
                }
                else{//Load physical positions
                    mapData->readMapDataVCF(vcfFilename, hapData->nloci, hapData->skipQueue);
                }
            }
            else
            {
                if (CALC_XP || CALC_XPNSL)
                {
                    hapData->initParams(UNPHASED, false, p.MAF, p.numThreads, flog, p.LOW_MEM);
                    hapData2->initParams(UNPHASED, false, p.MAF, p.numThreads, flog, p.LOW_MEM);

                    hapData->readHapData(hapFilename);
                    hapData2->readHapData(hapFilename2);

                    // handleData(hapFilename, "HAP", *hapData);
                    // handleData(hapFilename2, "HAP", *hapData2);

                    if (hapData->nloci != hapData2->nloci)
                    {
                        std::cerr << "ERROR: Haplotypes from " << hapFilename << " and " << hapFilename2 << " do not have the same number of loci.\n";
                        return 1;
                    }
                }else{
                    hapData->initParams(UNPHASED, p.SKIP, p.MAF, p.numThreads, flog, p.LOW_MEM);
                    hapData->readHapData(hapFilename);

                    // handleData(hapFilename, "HAP", *hapData);
                }
                mapData->readMapData(mapFilename, hapData->nloci, USE_PMAP, hapData->skipQueue);
            }
        }
        catch (exception& e)
        {
            cerr << "ERROR: " << e.what() << endl;
            return 1;
        }
        

        auto end_reading = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> read_duration =  end_reading - start_reading;
        //FLOG
        cout<<("Input file loaded in "+to_string(read_duration.count())+" s.")<<endl;

        (*(flog))<<("Input file loaded in "+to_string(read_duration.count())+" s.\n")<<endl;;
        //mapData.print();

        // Check if map is in order
        for (int i = 1; i < mapData->nloci; i++) {
            if ( mapData->mapEntries[i].physicalPos <  mapData->mapEntries[i-1].physicalPos ) {
                std::cerr << "ERROR: Variant physical position must be monotonically increasing.\n";
                std::cerr << "\t" << mapData->mapEntries[i].locusName << " " << mapData->mapEntries[i].physicalPos << " appears after";
                std::cerr << "\t" <<  mapData->mapEntries[i-1].locusName << " " << mapData->mapEntries[i-1].physicalPos << "\n";
                return 1;
            }
            if ( !CALC_NSL && mapData->mapEntries[i].geneticPos  < mapData->mapEntries[i-1].geneticPos  ) {
                std::cerr << "ERROR: Variant genetic position must be monotonically increasing.\n";
                std::cerr << "\t" << mapData->mapEntries[i].locusName << " " << mapData->mapEntries[i].geneticPos << " appears after";
                std::cerr << "\t" << mapData->mapEntries[i-1].locusName << " " << mapData->mapEntries[i-1].geneticPos << "\n";
                return 1;
            }
        }
        return 0;
    }

    /*
    bool is_gzipped(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        
        if (!file.is_open()) {
            std::cerr << "Error: Unable to open file " << filename << std::endl;
            return false;
        }

        // Read the first two bytes
        std::vector<unsigned char> buffer(2);
        file.read(reinterpret_cast<char*>(buffer.data()), buffer.size());

        // Close the file
        file.close();

        // Check if the first two bytes are the gzip magic numbers
        return buffer[0] == 0x1F && buffer[1] == 0x8B;
    }
    */


    /*
    void handleData(string filename, string EXT, HapData& hapData){
        if(EXT == "VCF"){
            //hapData.readHapDataVCF(filename);
            VCFSerialReader dataReader(filename, &hapData);
        }else if(EXT == "HAP"){
            //hapData.readHapData(filename);
            HapSerialReader dataReader(filename, &hapData);
        }else if(EXT == "TPED"){
            hapData.readHapDataTPED(filename);
        }

        if(false){
            if(EXT == "VCF"){
                hapData.readHapDataVCF(filename);
                // if(is_gzipped(filename)){
                //     cout<<"Gzipped file: so reading in serial"<<endl;
                //     VCFSerialReader dataReader(filename, &hapData);
                //     if(hapData.DEBUG_FLAG=="VCF") cout<<"Number of loci from serial reader: "<<hapData.nloci<<endl;
                // }else{
                //     cout<<"Plain text VCF file: so reading in parallel"<<endl;
                //     VCFParallelReader dataReader(filename, &hapData);
                // }
            }else if(EXT == "HAP"){
                //HapParallelReader dataReader(filename, &hapData);
                //HapSerialReader dataReader(filename, &hapData);
                hapData.readHapData(filename);

            }else if(EXT == "TPED"){
                hapData.readHapDataTPED(filename);
            }


            // // FINAL PHASE: FLIP if allowed
            if(hapData.benchmark_flag == "XOR" && hapData.benchmark_flag2 == "FLIP"){
                int count_flip = 0;
                for (int locus_after_filter = 0; locus_after_filter < hapData.nloci; locus_after_filter++){
                    if(hapData.hapEntries[locus_after_filter].positions.size() > hapData.nhaps/2){
                        hapData.hapEntries[locus_after_filter].flipped = true;
                        count_flip++;

                        vector<int> copy_pos;
                        copy_pos.reserve(hapData.nhaps - hapData.hapEntries[locus_after_filter].positions.size());
                        int cnt = 0;
                        for(int i = 0; i< hapData.nhaps; i++){
                            int curr =  hapData.hapEntries[locus_after_filter].positions[cnt];
                            if(i==curr){
                                cnt++;
                            }else{
                                copy_pos.push_back(i);
                            }
                        }
                        
                        hapData.hapEntries[locus_after_filter].positions = copy_pos;
                        vector<int>().swap(copy_pos);
                        // vector<int> zero_positions(this->nhaps - this->hapEntries[locus_after_filter].positions.size());
                        // int j = 0;
                        // int front_one = this->hapEntries[locus_after_filter].positions[j++];
                        // for(int i=0; i<nhaps; i++){
                        //     if(i==front_one){
                        //         front_one = this->hapEntries[locus_after_filter].positions[j++];
                        //     }else{
                        //         zero_positions.push_back(i);
                        //     }   
                        // }
                        // this->hapEntries[locus_after_filter].positions = zero_positions;
                    }else{
                        hapData.hapEntries[locus_after_filter].flipped = false;
                    }
                }
                cout<<"### Count of flipped loci: "<<count_flip<<endl;
            }

        }
        
        
    }
    */

};

#endif
