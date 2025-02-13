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

#include "selscan-cli.h"
#include "json.hpp"
#include <fstream>


void initalizeParameters(param_t &params,int argc, char *argv[]){
	params.setPreamble(PREAMBLE);
    params.addFlag(ARG_THREAD, DEFAULT_THREAD, "", HELP_THREAD);
    params.addFlag(ARG_FILENAME_POP1, DEFAULT_FILENAME_POP1, "", HELP_FILENAME_POP1);
    params.addFlag(ARG_FILENAME_POP2, DEFAULT_FILENAME_POP2, "", HELP_FILENAME_POP2);

    params.addFlag(ARG_FILENAME_POP1_THAP, DEFAULT_FILENAME_POP1_THAP, "", HELP_FILENAME_POP1_THAP);
    params.addFlag(ARG_FILENAME_POP2_THAP, DEFAULT_FILENAME_POP2_THAP, "", HELP_FILENAME_POP2_THAP);

    params.addFlag(ARG_FILENAME_POP1_TPED, DEFAULT_FILENAME_POP1_TPED, "", HELP_FILENAME_POP1_TPED);
    params.addFlag(ARG_FILENAME_POP2_TPED, DEFAULT_FILENAME_POP2_TPED, "", HELP_FILENAME_POP2_TPED);
    params.addFlag(ARG_FILENAME_POP1_VCF, DEFAULT_FILENAME_POP1_VCF, "", HELP_FILENAME_POP1_VCF);
    params.addFlag(ARG_FILENAME_POP2_VCF, DEFAULT_FILENAME_POP2_VCF, "", HELP_FILENAME_POP2_VCF);
    params.addFlag(ARG_FILENAME_MAP, DEFAULT_FILENAME_MAP, "", HELP_FILENAME_MAP);
    params.addFlag(ARG_PMAP, DEFAULT_PMAP, "", HELP_PMAP);
    params.addFlag(ARG_OUTFILE, DEFAULT_OUTFILE, "", HELP_OUTFILE);
    params.addFlag(ARG_CUTOFF, DEFAULT_CUTOFF, "", HELP_CUTOFF);
    params.addFlag(ARG_MAX_GAP, DEFAULT_MAX_GAP, "", HELP_MAX_GAP);
    params.addFlag(ARG_GAP_SCALE, DEFAULT_GAP_SCALE, "", HELP_GAP_SCALE);
    params.addFlag(ARG_IHS, DEFAULT_IHS, "", HELP_IHS);
    params.addFlag(ARG_XPNSL, DEFAULT_XPNSL, "", HELP_XPNSL);
    params.addFlag(ARG_UNPHASED, DEFAULT_UNPHASED, "", HELP_UNPHASED);
    params.addFlag(ARG_NSL, DEFAULT_NSL, "", HELP_NSL);
    params.addFlag(ARG_IHS_DETAILED, DEFAULT_IHS_DETAILED, "", HELP_IHS_DETAILED);
    params.addFlag(ARG_SOFT, DEFAULT_SOFT, "", HELP_SOFT);
    params.addFlag(ARG_XP, DEFAULT_XP, "", HELP_XP);
    params.addFlag(ARG_ALT, DEFAULT_ALT, "", HELP_ALT);
    params.addFlag(ARG_MAF, DEFAULT_MAF, "", HELP_MAF);
    params.addFlag(ARG_EHH, DEFAULT_EHH, "", HELP_EHH);

    params.addFlag(ARG_EHH12, DEFAULT_EHH12, "", HELP_EHH12);

    params.addFlag(ARG_QWIN, DEFAULT_QWIN, "", HELP_QWIN);
    //params.addFlag(ARG_SOFT_K, DEFAULT_SOFT_K, "SILENT", HELP_SOFT_K);
    params.addFlag(ARG_MAX_EXTEND, DEFAULT_MAX_EXTEND, "", HELP_MAX_EXTEND);
    params.addFlag(ARG_MAX_EXTEND_NSL, DEFAULT_MAX_EXTEND_NSL, "", HELP_MAX_EXTEND_NSL);
    params.addFlag(ARG_SKIP, DEFAULT_SKIP, "", HELP_SKIP);
    params.addFlag(ARG_KEEP, DEFAULT_KEEP, "", HELP_KEEP);
    params.addFlag(ARG_TRUNC, DEFAULT_TRUNC, "", HELP_TRUNC);
    params.addFlag(ARG_PI, DEFAULT_PI, "", HELP_PI);

    //params.addFlag(ARG_LOW_MEM, DEFAULT_LOW_MEM, "", HELP_LOW_MEM);
    params.addFlag(ARG_MISSING_FLAG, DEFAULT_MISSING_FLAG, "", HELP_MISSING_FLAG);
    params.addFlag(ARG_IMPUTE_FLAG, DEFAULT_IMPUTE_FLAG, "", HELP_IMPUTE_FLAG);

    //params.addFlag(ARG_MULTI_CHR, DEFAULT_MULTI_CHR, "", HELP_MULTI_CHR);

    params.addFlag(ARG_PI_WIN, DEFAULT_PI_WIN, "", HELP_PI_WIN);
    params.addFlag(ARG_WAGH, DEFAULT_WAGH, "", HELP_WAGH);

    params.addFlag(ARG_MULTI_PARAMS, DEFAULT_MULTI_PARAMS, "", HELP_MULTI_PARAMS);


    try
    {
        params.parseCommandLine(argc, argv);
    }
    catch (...)
    {
        throw runtime_error("Could not parse command line arguments.\n");
    }
}


void getBaseParamFromCmdLine(param_t& params, param_main &p){
    p.hapFilename = params.getStringFlag(ARG_FILENAME_POP1);
    p.hapFilename2 = params.getStringFlag(ARG_FILENAME_POP2);
    p.mapFilename = params.getStringFlag(ARG_FILENAME_MAP);
    p.tpedFilename = params.getStringFlag(ARG_FILENAME_POP1_TPED);
    p.tpedFilename2 = params.getStringFlag(ARG_FILENAME_POP2_TPED);
    p.vcfFilename = params.getStringFlag(ARG_FILENAME_POP1_VCF);
    p.vcfFilename2 = params.getStringFlag(ARG_FILENAME_POP2_VCF);
    
    p.thapFilename = params.getStringFlag(ARG_FILENAME_POP1_THAP);
    p.thapFilename2 = params.getStringFlag(ARG_FILENAME_POP2_THAP);
    
    p.jsonFilename = params.getStringFlag(ARG_MULTI_PARAMS);

    p.TPED = false;
    if (p.tpedFilename.compare(DEFAULT_FILENAME_POP1_TPED) != 0) p.TPED = true;

    p.VCF = false;
    if (p.vcfFilename.compare(DEFAULT_FILENAME_POP1_VCF) != 0) p.VCF = true;

    p.THAP = false;
    if (p.hapFilename.compare(DEFAULT_FILENAME_POP1_THAP) != 0) p.THAP = true;


    p.outFilename = params.getStringFlag(ARG_OUTFILE);
    p.query = params.getStringFlag(ARG_EHH);
    p.query_ehh12 = params.getStringFlag(ARG_EHH12);

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
    p.SINGLE_EHH12 = false;

    // p.EHH1K = params.getIntFlag(ARG_SOFT_K);

    int queryLoc = -1;
    if (p.query.compare(DEFAULT_EHH) != 0) p.SINGLE_EHH = true;
    if (p.query_ehh12.compare(DEFAULT_EHH12) != 0) p.SINGLE_EHH12 = true;


    p.QWIN = params.getIntFlag(ARG_QWIN);

    p.CALC_PI = params.getBoolFlag(ARG_PI);
    p.PI_WIN = params.getIntFlag(ARG_PI_WIN);

    //params.getBoolFlag(ARG_LOW_MEM);

    p.MISSING_ALLOWED = params.getBoolFlag(ARG_MISSING_FLAG);

    p.MAX_EXTEND_NSL = params.getIntFlag(ARG_MAX_EXTEND_NSL); //p.MAX_EXTEND_NSL = ( params.getIntFlag(ARG_MAX_EXTEND_NSL) <= 0 ) ? hm.mapData.nloci : params.getIntFlag(ARG_MAX_EXTEND_NSL);
    p.MAX_EXTEND =  params.getIntFlag(ARG_MAX_EXTEND);
    p.TRUNC = params.getBoolFlag(ARG_TRUNC);


    // p.CHR_LIST = params.getStringFlag(ARG_MULTI_CHR);
    // p.MULTI_CHR = false;
    // if (p.CHR_LIST.compare(DEFAULT_MULTI_CHR) != 0) p.MULTI_CHR = true;

 
    if(params.getBoolFlag(ARG_IMPUTE_FLAG)){
        p.MISSING_MODE = "NO_IMPUTE";
    }else{
        p.MISSING_MODE = "ONE_IMPUTE";
    }

    p.MULTI_MAF = false;
    p.MULTI_PARAMS = false;
    if (p.jsonFilename.compare(DEFAULT_MULTI_PARAMS) != 0) {
        p.MULTI_PARAMS = true;
    }

    p.SKIP = !params.getBoolFlag(ARG_KEEP);//params.getBoolFlag(ARG_SKIP);
    
}

//Returns 1 if error
//bool checkParameters(param_t &params,int argc, char *argv[]){
void checkParameters(param_main &p){

    // string hapFilename = p.hapFilename;
    // string hapFilename2 = p.hapFilename2;
    // string mapFilename = p.mapFilename;
    // string tpedFilename = p.tpedFilename;
    // string tpedFilename2 = p.tpedFilename2;
    // string vcfFilename = p.vcfFilename;
    // string vcfFilename2 = p.vcfFilename2;
    // bool TPED = p.TPED;
    // bool VCF = p.VCF;

    // bool TPED = false;
    // if (tpedFilename.compare(DEFAULT_FILENAME_POP1_TPED) != 0) TPED = true;

    // bool VCF = false;
    // if (vcfFilename.compare(DEFAULT_FILENAME_POP1_VCF) != 0) VCF = true;

    if (p.VCF && p.TPED) {
        throw runtime_error("Please choose only one of TPED, VCF, or HAP formatted files.\n");
    }

    if ( (p.VCF || p.TPED) && (p.hapFilename.compare(DEFAULT_FILENAME_POP1) != 0 || p.hapFilename2.compare(DEFAULT_FILENAME_POP2) != 0) ) {
        throw runtime_error("Please choose only one of TPED, VCF, or HAP formatted files.\n");
    }
    
	// string outFilename = params.getStringFlag(ARG_OUTFILE);
    // string query = params.getStringFlag(ARG_EHH);


    //bool TRUNC = params.getBoolFlag(ARG_TRUNC);

    // int EHH1K = params.getIntFlag(ARG_SOFT_K);


    // bool TRUNC = p.TRUNC;
    // int EHH1K = 2; //params.getIntFlag(ARG_SOFT_K);
    // int QWIN = p.QWIN;

    // bool CALC_PI = p.CALC_PI;
    // int PI_WIN = p.PI_WIN;

    //TODON
    char PI_WIN_str[50];
    snprintf(PI_WIN_str, 50, "%d", p.PI_WIN);


    //DEBUG::: cout<<  CALC_XPNSL << " " << CALC_IHS << " " << CALC_XP << " " << SINGLE_EHH << " " << CALC_PI << " " << CALC_NSL << " " << CALC_SOFT << endl;
    
    //     if (CALC_XPNSL + CALC_IHS + CALC_XP + SINGLE_EHH + CALC_PI + CALC_NSL + CALC_SOFT != 1)
    // {
    //     cerr << "ERROR: Must specify one and only one of \n\tEHH (" << ARG_EHH
    //          << ")\n\tiHS (" << ARG_IHS
    //          << ")\n\tXP-EHH (" << ARG_XP
    //          << ")\n\tPI (" << ARG_PI
    //          << ")\n\tnSL (" << ARG_NSL
    //          << ")\n\tXP-nSL (" << ARG_XPNSL
    //          << ")\n\tiHH12 (" << ARG_SOFT
    //          << ")\n";
    //     return 1;
    // }

        if (p.CALC_XPNSL + p.CALC_IHS + p.CALC_XP + p.SINGLE_EHH + p.CALC_PI + p.CALC_NSL + p.CALC_SOFT + p.SINGLE_EHH12 < 1)
    {
        // cerr << "ERROR: Must specify one of \n\tEHH (" << ARG_EHH
        //      << ")\n\tiHS (" << ARG_IHS
        //      << ")\n\tXP-EHH (" << ARG_XP
        //      << ")\n\tPI (" << ARG_PI
        //      << ")\n\tnSL (" << ARG_NSL
        //      << ")\n\tXP-nSL (" << ARG_XPNSL
        //      << ")\n\tiHH12 (" << ARG_SOFT
        //      << ")\n";
        throw runtime_error("Must specify one of \n\tEHH (" + ARG_EHH + ")\n\tiHS (" + ARG_IHS + ")\n\tXP-EHH (" + ARG_XP + ")\n\tPI (" + ARG_PI + ")\n\tnSL (" + ARG_NSL + ")\n\tXP-nSL (" + ARG_XPNSL + ")\n\tiHH12 (" + ARG_SOFT + ")\n\tEHH12 (" + ARG_EHH12 + ")\n");
    }

    if (p.CALC_XPNSL + p.CALC_IHS + p.CALC_XP + p.SINGLE_EHH + p.CALC_PI + p.CALC_NSL + p.CALC_SOFT + p.SINGLE_EHH12> 1)
    {
        cerr<<"WARNING: Running multiple statistics at once. Make sure parameters are consistent across runs.\n";
        p.MULTI_MAF = true;
    }

    if (p.WRITE_DETAILED_IHS && !p.CALC_IHS) {
        throw runtime_error("The flag " + ARG_IHS_DETAILED + " must be used with " + ARG_IHS + " \n");
    }

    /*
        if (SINGLE_EHH && CALC_XP)
        {
            cerr << "Single query with XP-EHH is not yet available.\n";
            return 1;
        }
    */

    if(p.MISSING_ALLOWED){
        throw runtime_error("--missing flag is still under development, so it is disabled in this version.\n");
    }

    if (p.WAGH && p.UNPHASED){
        string error_msg = "--wagh and --unphased currently incompatible.\n\t\
        Consider --xpehh or --xpnsl with --unphased for two population selection statistics.\n";
        throw runtime_error(error_msg);
    }

    if (p.CALC_PI && p.UNPHASED){
        throw runtime_error("--pi and --unphased currently incompatible.\n");
    }

    if (p.CALC_SOFT && p.UNPHASED){
        throw runtime_error("--ihh12 and --unphased currently incompatible.\n");
    }
    if (p.numThreads < 1)
    {
        throw runtime_error("Must run with one or more threads.\n");
    }
    if (p.SCALE_PARAMETER < 1)
    {
        throw runtime_error("Scale parameter must be positive.\n");
    }
    if (p.MAX_GAP < 1)
    {
        throw runtime_error("Max gap parameter must be positive.\n");
    }
    if (p.EHH_CUTOFF <= 0 || p.EHH_CUTOFF >= 1)
    {
        throw runtime_error("EHH cut off must be > 0 and < 1.\n");
    }
    if (p.TPED)
    {
        if ((!p.CALC_XP && !p.CALC_XPNSL) && p.tpedFilename2.compare(DEFAULT_FILENAME_POP2_TPED) != 0)
        {
            throw runtime_error("You are not calculating XP stats but have given a second data file (" + p.tpedFilename2 + ").\n");
        }
    }
    else if (p.VCF) {
        if ((!p.CALC_XP && !p.CALC_XPNSL) && p.vcfFilename2.compare(DEFAULT_FILENAME_POP2_VCF) != 0)
        {
            throw runtime_error("You are not calculating XP stats but have given a second data file (" + p.vcfFilename2 + ").\n");
        }

        if ((!p.CALC_NSL && !p.CALC_XPNSL) && (p.mapFilename.compare(DEFAULT_FILENAME_MAP) == 0 && !p.USE_PMAP)) {
            throw runtime_error("Must also provide a mapfile.\n");
        }
    }
    else if(p.THAP)
    {
        if ((!p.CALC_XP && !p.CALC_XPNSL) && p.thapFilename2.compare(DEFAULT_FILENAME_POP2_THAP) != 0)
        {
            throw runtime_error("You are not calculating XP stats but have given a second data file (" + p.thapFilename2 + ").\n");
        }
        if (p.mapFilename.compare(DEFAULT_FILENAME_MAP) == 0) {
            throw runtime_error("Must also provide a mapfile.\n");
        }
    }else{
        if ((!p.CALC_XP && !p.CALC_XPNSL) && p.hapFilename2.compare(DEFAULT_FILENAME_POP2) != 0)
        {
            throw runtime_error("You are not calculating XP stats but have given a second data file (" + p.hapFilename2 + ").\n");
        }
        if (p.mapFilename.compare(DEFAULT_FILENAME_MAP) == 0) {
            throw runtime_error("Must also provide a mapfile.\n");
        }
    }
    // if (EHH1K < 1)
    // {
    //     cerr << "ERROR: EHH1K must be > 0.\n";
    //     return 1;
    // }

    if (p.PI_WIN < 1)
    {
        throw runtime_error("pi window must be > 0.\n");
    }

    // if (p.MULTI_CHR && !p.VCF)
    // {
    //     throw runtime_error("--multi-chr flag only works with VCF files.\n");
    // }

    //check that string is comma-separated list
    // string CHR_LIST = p.CHR_LIST;
    // if (CHR_LIST.length()==0)
    // {
    //     cerr << "ERROR: you must provide a comma-separated list after --chr flag. \n";
    //     return 1;
    // }

    // if (p.MULTI_PARAMS && p.jsonFilename.compare(DEFAULT_MULTI_PARAMS) == 0) {
    //     throw runtime_error("Must provide a JSON file.\n");
    // }
    
    if(p.SKIP){
        cerr << "WARNING: " << ARG_SKIP << " is now on by dafault.  This flag no longer has a function.\n";
    }
    // if(params.getBoolFlag(ARG_SKIP)){
    //     cerr << "WARNING: " << ARG_SKIP << " is now on by dafault.  This flag no longer has a function.\n";
    // }
    //return 0;
}

void jsonParse(param_main &base_p, vector<param_main> &multi_params){
    std::ifstream inputFile(base_p.jsonFilename);

    if (!inputFile.is_open()) {
        cerr<<"ERROR: Failed to open json file "<< base_p.jsonFilename << " containing parameters.\n";
        exit(1);
    }

    nlohmann::json data;
    inputFile >> data;
    inputFile.close();
    
    // Iterate over the array of params
    for (const auto& param : data["params"]) {
        param_main p(base_p); // make a copy of the base parameters
        //cout<<ARG_MAX_EXTEND<<" "<<ARG_MAF<<" "<<ARG_MAX_EXTEND_NSL<<endl;
        string MAX_EXTEND = ARG_MAX_EXTEND.substr(2);
        string MAF = ARG_MAF.substr(2);
        string MAX_EXTEND_NSL = ARG_MAX_EXTEND_NSL.substr(2);
        string CUTOFF = ARG_CUTOFF.substr(2);
        // all other params are ignored

        string FILENAME_POP1 = ARG_FILENAME_POP1.substr(2);
        string FILENAME_POP2 = ARG_FILENAME_POP2.substr(2);

        if (param.contains(MAF)) { // Check if the field exists and retrieve it if it does
            p.MAF = param[MAF];
            cerr<<"DEBUG:::"<<ARG_MAF<<" "<<p.MAF<<endl;
        }

        if (param.contains(MAX_EXTEND)) { // Check if the field exists and retrieve it if it does
            p.MAX_EXTEND = param[MAX_EXTEND];
            cerr<<"DEBUG:::"<<ARG_MAX_EXTEND<<" "<<p.MAX_EXTEND<<endl;
        }

        if (param.contains(MAX_EXTEND_NSL)) { // Check if the field exists and retrieve it if it does
            p.MAX_EXTEND_NSL = param[MAX_EXTEND_NSL];
            cerr<<"DEBUG:::"<<ARG_MAX_EXTEND_NSL<<" "<<p.MAX_EXTEND_NSL<<endl;
        }
        if(param.contains(CUTOFF)){
            p.EHH_CUTOFF = param[CUTOFF];
            cerr<<"DEBUG:::"<<ARG_CUTOFF<<" "<<p.EHH_CUTOFF<<endl;
        }
        checkParameters(p);
        multi_params.push_back(p);
    }
}