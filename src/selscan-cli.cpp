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

//Returns 1 if error
bool initalizeParameters(param_t &params,int argc, char *argv[]){
	params.setPreamble(PREAMBLE);
    params.addFlag(ARG_THREAD, DEFAULT_THREAD, "", HELP_THREAD);
    params.addFlag(ARG_FILENAME_POP1, DEFAULT_FILENAME_POP1, "", HELP_FILENAME_POP1);
    params.addFlag(ARG_FILENAME_POP2, DEFAULT_FILENAME_POP2, "", HELP_FILENAME_POP2);
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
    params.addFlag(ARG_MULTICHR_FLAG, DEFAULT_MULTICHR_FLAG, "", HELP_MULTICHR_FLAG);
    params.addFlag(ARG_CHR_LIST, DEFAULT_CHR_LIST, "", HELP_CHR_LIST);


    params.addFlag(ARG_PI_WIN, DEFAULT_PI_WIN, "", HELP_PI_WIN);
    params.addFlag(ARG_WAGH, DEFAULT_WAGH, "", HELP_WAGH);

    try
    {
        params.parseCommandLine(argc, argv);
    }
    catch (...)
    {
        return 1;
    }

    return 0;
}

//Returns 1 if error
bool checkParameters(param_t &params,int argc, char *argv[]){

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

    if (VCF && TPED) {
        cerr << "ERROR: Please choose only one of TPED, VCF, or HAP formatted files.\n";
        return 1;
    }

    if ( (VCF || TPED) && (hapFilename.compare(DEFAULT_FILENAME_POP1) != 0 || hapFilename2.compare(DEFAULT_FILENAME_POP2) != 0) ) {
        cerr << "ERROR: Please choose only one of TPED, VCF, or HAP formatted files.\n";
        return 1;
    }

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
    bool SINGLE_EHH =  false;
    if (query.compare(DEFAULT_EHH) != 0) SINGLE_EHH = true;

    bool SKIP = !params.getBoolFlag(ARG_KEEP);//params.getBoolFlag(ARG_SKIP);
    if(params.getBoolFlag(ARG_SKIP)){
        cerr << "WARNING: " << ARG_SKIP << " is now on by dafault.  This flag no longer has a function.\n";
    }
    //bool TRUNC = params.getBoolFlag(ARG_TRUNC);

    // int EHH1K = params.getIntFlag(ARG_SOFT_K);

    bool CALC_PI = params.getBoolFlag(ARG_PI);
    int PI_WIN = params.getIntFlag(ARG_PI_WIN);


    cout<<  CALC_XPNSL << " " << CALC_IHS << " " << CALC_XP << " " << SINGLE_EHH << " " << CALC_PI << " " << CALC_NSL << " " << CALC_SOFT << endl;
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

    if (WRITE_DETAILED_IHS && !CALC_IHS) {
        cerr << "ERROR: The flag " << ARG_IHS_DETAILED << " must be used with " << ARG_IHS << " \n";
        return 1;
    }

    /*
        if (SINGLE_EHH && CALC_XP)
        {
            cerr << "Single query with XP-EHH is not yet available.\n";
            return 1;
        }
    */

    if (WAGH && UNPHASED){
        cerr << "ERROR: --wagh and --unphased currently incompatible.\n\t\
        Consider --xpehh or --xpnsl with --unphased for two population selection statistics.\n";
        return 1;
    }

    if (CALC_PI && UNPHASED){
        cerr << "ERROR: --pi and --unphased currently incompatible.\n";
        return 1;
    }

    if (CALC_SOFT && UNPHASED){
        cerr << "ERROR: --ihh12 and --unphased currently incompatible.\n";
        return 1;
    }
    if (numThreads < 1)
    {
        cerr << "ERROR: Must run with one or more threads.\n";
        return 1;
    }
    if (SCALE_PARAMETER < 1)
    {
        cerr << "ERROR: Scale parameter must be positive.\n";
        return 1;
    }
    if (MAX_GAP < 1)
    {
        cerr << "ERROR: Max gap parameter must be positive.\n";
        return 1;
    }
    if (EHH_CUTOFF <= 0 || EHH_CUTOFF >= 1)
    {
        cerr << "ERROR: EHH cut off must be > 0 and < 1.\n";
        return 1;
    }
    if (TPED)
    {
        if ((!CALC_XP && !CALC_XPNSL) && tpedFilename2.compare(DEFAULT_FILENAME_POP2_TPED) != 0)
        {
            cerr << "ERROR: You are not calculating XP stats but have given a second data file (" << tpedFilename2 << ").\n";
            return 1;
        }
    }
    else if (VCF) {
        if ((!CALC_XP && !CALC_XPNSL) && vcfFilename2.compare(DEFAULT_FILENAME_POP2_VCF) != 0)
        {
            cerr << "ERROR: You are not calculating XP stats but have given a second data file (" << vcfFilename2 << ").\n";
            return 1;
        }

        if ((!CALC_NSL && !CALC_XPNSL) && (mapFilename.compare(DEFAULT_FILENAME_MAP) == 0 && !USE_PMAP)) {
            cerr << "ERROR: Must also provide a mapfile.\n";
            return 1;
        }
    }
    else
    {
        if ((!CALC_XP && !CALC_XPNSL) && hapFilename2.compare(DEFAULT_FILENAME_POP2) != 0)
        {
            cerr << "ERROR: You are not calculating XP stats but have given a second data file (" << hapFilename2 << ").\n";
            return 1;
        }
        if (mapFilename.compare(DEFAULT_FILENAME_MAP) == 0) {
            cerr << "ERROR: Must also provide a mapfile.\n";
            return 1;
        }

    }
    // if (EHH1K < 1)
    // {
    //     cerr << "ERROR: EHH1K must be > 0.\n";
    //     return 1;
    // }

    if (PI_WIN < 1)
    {
        cerr << "ERROR: pi window must be > 0.\n";
        return 1;
    }



    bool MULTICHR = params.getBoolFlag(ARG_MULTICHR_FLAG);
    if (MULTICHR && !VCF)
    {
        cerr << "ERROR: --multi-chr flag only works with VCF files.\n";
        return 1;
    }

    //check that string is comma-separated list
    string chrList = params.getStringFlag(ARG_CHR_LIST);
    if (chrList.length()==0)
    {
        cerr << "ERROR: you must provide a comma-separated list after --chr flag. \n";
        return 1;
    }
    
    return 0;
}



