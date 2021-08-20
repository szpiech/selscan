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
#include <pthread.h>
#include <map>
#include "selscan-data.h"
#include "selscan-pbar.h"
#include "binom.h"
#include "param_t.h"
#include "hamming_t.h"

using namespace std;

const string VERSION = "1.4.0";

const string PREAMBLE = "\nselscan v" + VERSION + " -- a program to calculate EHH-based scans for positive selection in genomes.\n\
Source code and binaries can be found at <https://www.github.com/szpiech/selscan>.\n\
\n\
selscan currently implements EHH, iHS, XP-EHH, and nSL.\n\
\n\
Citations:\n\
\n\
selscan: ZA Szpiech and RD Hernandez (2014) MBE 31: 2824-2827.\n\
iHH12: R Torres et al. (2018) PLoS Genetics 15: e1007898.\n\
       N Garud et al. (2015) PLoS Genetics 11: 1–32.\n\
nSL: A Ferrer-Admetlla et al. (2014) MBE 31: 1275-1291.\n\
XP-nSL: Szpiech et al. (2021) Evol Lett 5: 408-421.\n\
XP-EHH: PC Sabeti et al. (2007) Nature 449: 913–918.\n\
        K Wagh et al. (2012) PloS ONE 7: e44751.\n\
iHS: BF Voight et al. (2006) PLoS Biology 4: e72.\n\
EHH: PC Sabeti et al. (2002) Nature 419: 832–837.\n\
\n\
To calculate EHH:\n\
\n\
./selscan --ehh <locusID> --vcf <vcf> --map <mapfile> --out <outfile>\n\
\n\
To calculate iHS:\n\
\n\
./selscan --ihs --vcf <vcf> --map <mapfile> --out <outfile>\n\
\n\
To calculate nSL:\n\
\n\
./selscan --nsl --vcf <vcf> --out <outfile>\n\
\n\
To calculate XP-nSL:\n\
\n\
./selscan --xpnsl --vcf <vcf> --vcf-ref <vcf> --out <outfile>\n\
\n\
To calculate iHH12:\n\
\n\
./selscan --ihh12 --vcf <vcf> --map <mapfile> --out <outfile>\n\
\n\
To calculate XP-EHH:\n\
\n\
./selscan --xpehh --vcf <vcf> --vcf-ref <vcf> --map <mapfile> --out <outfile>\n";

const string ARG_THREAD = "--threads";
const int DEFAULT_THREAD = 1;
const string HELP_THREAD = "The number of threads to spawn during the calculation.\n\
\tPartitions loci across threads.";

const string ARG_FILENAME_POP1_TPED = "--tped";
const string DEFAULT_FILENAME_POP1_TPED = "__hapfile1";
const string HELP_FILENAME_POP1_TPED = "A TPED file containing haplotype and map data.\n\
\tVariants should be coded 0/1";

const string ARG_FILENAME_POP2_TPED = "--tped-ref";
const string DEFAULT_FILENAME_POP2_TPED = "__hapfile2";
const string HELP_FILENAME_POP2_TPED = "A TPED file containing haplotype and map data.\n\
\tVariants should be coded 0/1. This is the 'reference'\n\
\tpopulation for XP-EHH calculations and should contain the same number\n\
\tof loci as the query population. Ignored otherwise.";

const string ARG_FILENAME_POP1_VCF = "--vcf";
const string DEFAULT_FILENAME_POP1_VCF = "__hapfile1";
const string HELP_FILENAME_POP1_VCF = "A VCF file containing haplotype data.\n\
\tA map file must be specified with --map.";

const string ARG_FILENAME_POP2_VCF = "--vcf-ref";
const string DEFAULT_FILENAME_POP2_VCF = "__hapfile2";
const string HELP_FILENAME_POP2_VCF = "A VCF file containing haplotype and map data.\n\
\tVariants should be coded 0/1. This is the 'reference'\n\
\tpopulation for XP-EHH calculations and should contain the same number\n\
\tof loci as the query population. Ignored otherwise.";

const string ARG_FILENAME_POP1 = "--hap";
const string DEFAULT_FILENAME_POP1 = "__hapfile1";
const string HELP_FILENAME_POP1 = "A hapfile with one row per haplotype, and one column per\n\
\tvariant. Variants should be coded 0/1";

const string ARG_FILENAME_POP2 = "--ref";
const string DEFAULT_FILENAME_POP2 = "__hapfile2";
const string HELP_FILENAME_POP2 = "A hapfile with one row per haplotype, and one column per\n\
\tvariant. Variants should be coded 0/1. This is the 'reference'\n\
\tpopulation for XP-EHH calculations.  Ignored otherwise.";

const string ARG_FILENAME_MAP = "--map";
const string DEFAULT_FILENAME_MAP = "__mapfile";
const string HELP_FILENAME_MAP = "A mapfile with one row per variant site.\n\
\tFormatted <chr#> <locusID> <genetic pos> <physical pos>.";

const string ARG_PMAP = "--pmap";
const bool DEFAULT_PMAP = false;
const string HELP_PMAP = "Use physical map instead of a genetic map.";

const string ARG_OUTFILE = "--out";
const string DEFAULT_OUTFILE = "outfile";
const string HELP_OUTFILE = "The basename for all output files.";

const string ARG_CUTOFF = "--cutoff";
const double DEFAULT_CUTOFF = 0.05;
const string HELP_CUTOFF = "The EHH decay cutoff.";

const string ARG_MAX_GAP = "--max-gap";
const int DEFAULT_MAX_GAP = 200000;
const string HELP_MAX_GAP = "Maximum allowed gap in bp between two snps.";

const string ARG_GAP_SCALE = "--gap-scale";
const int DEFAULT_GAP_SCALE = 20000;
const string HELP_GAP_SCALE = "Gap scale parameter in bp. If a gap is encountered between\n\
\ttwo snps > GAP_SCALE and < MAX_GAP, then the genetic distance is\n\
\tscaled by GAP_SCALE/GAP.";

const string ARG_IHS = "--ihs";
const bool DEFAULT_IHS = false;
const string HELP_IHS = "Set this flag to calculate iHS.";

const string ARG_IHS_DETAILED = "--ihs-detail";
const bool DEFAULT_IHS_DETAILED = false;
const string HELP_IHS_DETAILED = "Set this flag to write out left and right iHH scores for '1' and '0' in addition to iHS.";

const string ARG_UNPHASED = "--unphased";
const bool DEFAULT_UNPHASED = false;
const string HELP_UNPHASED = "Set this flag to use multilocus genotypes.";

const string ARG_XPNSL = "--xpnsl";
const bool DEFAULT_XPNSL = false;
const string HELP_XPNSL = "Set this flag to calculate XP-nSL.";

const string ARG_NSL = "--nsl";
const bool DEFAULT_NSL = false;
const string HELP_NSL = "Set this flag to calculate nSL.";

const string ARG_SOFT = "--ihh12";
const bool DEFAULT_SOFT = false;
const string HELP_SOFT = "Set this flag to calculate iHH12.";

const string ARG_SOFT_K = "--k";
const int DEFAULT_SOFT_K = 2;
const string HELP_SOFT_K = "Specify K to compute for iHH1K.";

const string ARG_XP = "--xpehh";
const bool DEFAULT_XP = false;
const string HELP_XP = "Set this flag to calculate XP-EHH.";

const string ARG_ALT = "--alt";
const bool DEFAULT_ALT = false;
const string HELP_ALT = "Set this flag to calculate homozygosity based on the sum of the\n\
\tsquared haplotype frequencies in the observed data instead of using\n\
\tbinomial coefficients.";

const string ARG_WAGH = "--wagh";
const bool DEFAULT_WAGH = false;
const string HELP_WAGH = "Set this flag to calculate XP-EHH using definition of EHH which\n\
\tseparates core SNP alleles in the denominator.";

const string ARG_MAF = "--maf";
const double DEFAULT_MAF = 0.05;
const string HELP_MAF = "If a site has a MAF below this value, the program will not use\n\
\tit as a core snp.";

const string ARG_EHH = "--ehh";
const string DEFAULT_EHH = "__NO_LOCUS__";
const string HELP_EHH = "Calculate EHH of the '1' and '0' haplotypes at the specified\n\
\tlocus. Output: <physical dist> <genetic dist> <'1' EHH> <'0' EHH>";

const string ARG_QWIN = "--ehh-win";
const int DEFAULT_QWIN = 100000;
const string HELP_QWIN = "When calculating EHH, this is the length of the window in bp\n\
\tin each direction from the query locus.";

const string ARG_MAX_EXTEND = "--max-extend";
const int DEFAULT_MAX_EXTEND = 1000000;
const string HELP_MAX_EXTEND = "The maximum distance an EHH decay curve is allowed to extend from the core.\n\
\tSet <= 0 for no restriction.";

const string ARG_MAX_EXTEND_NSL = "--max-extend-nsl";
const int DEFAULT_MAX_EXTEND_NSL = 100;
const string HELP_MAX_EXTEND_NSL = "The maximum distance an nSL haplotype is allowed to extend from the core.\n\
\tSet <= 0 for no restriction.";

const string ARG_SKIP = "--skip-low-freq";
const bool DEFAULT_SKIP = false;
const string HELP_SKIP = "**This flag is now on by default. If you want to include low frequency variants\n\
in the construction of your haplotypes please use the --keep-low-freq flag.";

const string ARG_KEEP = "--keep-low-freq";
const bool DEFAULT_KEEP = false;
const string HELP_KEEP = "Include low frequency variants in the construction of your haplotypes.";

const string ARG_TRUNC = "--trunc-ok";
const bool DEFAULT_TRUNC = false;
const string HELP_TRUNC = "If an EHH decay reaches the end of a sequence before reaching the cutoff,\n\
\tintegrate the curve anyway (iHS and XPEHH only).\n\
\tNormal function is to disregard the score for that core.";

const string ARG_PI = "--pi";
const bool DEFAULT_PI = false;
const string HELP_PI = "Set this flag to calculate mean pairwise sequence difference in a sliding window.";

const string ARG_PI_WIN = "--pi-win";
const int DEFAULT_PI_WIN = 100;
const string HELP_PI_WIN = "Sliding window size in bp for calculating pi.";

pthread_mutex_t mutex_log = PTHREAD_MUTEX_INITIALIZER;

struct work_order_t
{
    int queryLoc;
    int id;

    string filename;

    HaplotypeData *hapData;

    HaplotypeData *hapData1;
    HaplotypeData *hapData2;

    MapData *mapData;

    double (*calc)(map<string, int> &, int, bool);

    double *ihs;
    double *ihhDerivedLeft;
    double *ihhDerivedRight;
    double *ihhAncestralLeft;
    double *ihhAncestralRight;
    double *freq;

    double *ihh1;
    double *freq1;

    double *ihh2;
    double *freq2;

    ofstream *flog;
    ofstream *fout;
    Bar *bar;

    param_t *params;
};

struct triplet_t
{
    double h1;
    double h12;
    double h2dh1;
};

void calculatePi(HaplotypeData *hapData, MapData *mapData, int winsize, string outFilename);

triplet_t calculateSoft(map<string, int> &count, int total);

void query_locus(void *work_order);
void query_locus_soft(void *order);

void calc_ihs(void *work_order);
void calc_nsl(void *work_order);
void calc_xpihh(void *work_order);
void calc_soft_ihs(void *order);

double calcFreq(HaplotypeData *hapData, int locus, bool unphased);
int queryFound(MapData *mapData, string query);
void fillColors(int **hapColor, map<string, int> &hapCount,
                string *haplotypeList, int hapListLength,
                int currentLoc, int &currentColor, bool left);
bool familyDidSplit(const string &hapStr, const int hapCount,
                    int **hapColor, const int nhaps, const int colorIndex,
                    const int previousLoc, string &mostCommonHap);

double calculateHomozygosity_Wagh(map<string, int> &count, int total, int derivedCount);

double calculateHomozygosity(map<string, int> &count, int total, bool ALT);

int main(int argc, char *argv[])
{
    cerr << "selscan v" + VERSION + "\n";
#ifdef PTW32_STATIC_LIB
    pthread_win32_process_attach_np();
#endif

    param_t params;
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
    params.addFlag(ARG_SOFT_K, DEFAULT_SOFT_K, "SILENT", HELP_SOFT_K);
    params.addFlag(ARG_MAX_EXTEND, DEFAULT_MAX_EXTEND, "", HELP_MAX_EXTEND);
    params.addFlag(ARG_MAX_EXTEND_NSL, DEFAULT_MAX_EXTEND_NSL, "", HELP_MAX_EXTEND_NSL);
    params.addFlag(ARG_SKIP, DEFAULT_SKIP, "", HELP_SKIP);
    params.addFlag(ARG_KEEP, DEFAULT_KEEP, "", HELP_KEEP);
    params.addFlag(ARG_TRUNC, DEFAULT_TRUNC, "", HELP_TRUNC);
    params.addFlag(ARG_PI, DEFAULT_PI, "", HELP_PI);
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
    sprintf(PI_WIN_str, "%d", PI_WIN);

    if (query.compare(DEFAULT_EHH) != 0) SINGLE_EHH = true;


    if (CALC_XPNSL + CALC_IHS + CALC_XP + SINGLE_EHH + CALC_PI + CALC_NSL + CALC_SOFT != 1)
    {
        cerr << "ERROR: Must specify one and only one of \n\tEHH (" << ARG_EHH
             << ")\n\tiHS (" << ARG_IHS
             << ")\n\tXP-EHH (" << ARG_XP
             << ")\n\tPI (" << ARG_PI
             << ")\n\tnSL (" << ARG_NSL
             << ")\n\tXP-nSL (" << ARG_XPNSL
             << ")\n\tiHH12 (" << ARG_SOFT
             << ")\n";
        return 1;
    }

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

    if (SINGLE_EHH) outFilename += ".ehh." + query;
    else if (CALC_IHS) outFilename += ".ihs";
    else if (CALC_NSL) outFilename += ".nsl";
    else if (CALC_XPNSL) outFilename += ".xpnsl";
    else if (CALC_XP) outFilename += ".xpehh";
    else if (CALC_SOFT) outFilename += ".ihh12";
    else if (CALC_PI) outFilename += ".pi." + string(PI_WIN_str) + "bp";

    if (ALT) outFilename += ".alt";

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
    if (EHH1K < 1)
    {
        cerr << "ERROR: EHH1K must be > 0.\n";
        return 1;
    }

    if (PI_WIN < 1)
    {
        cerr << "ERROR: pi window must be > 0.\n";
        return 1;
    }

    HaplotypeData *hapData, *hapData2;
    MapData *mapData;

    try
    {
        if (TPED)
        {
            hapData = readHaplotypeDataTPED(tpedFilename,UNPHASED);
            if (CALC_XP || CALC_XPNSL)
            {
                hapData2 = readHaplotypeDataTPED(tpedFilename2,UNPHASED);
                if (hapData->nloci != hapData2->nloci)
                {
                    cerr << "ERROR: Haplotypes from " << tpedFilename << " and " << tpedFilename2 << " do not have the same number of loci.\n";
                    return 1;
                }
            }
            mapData = readMapDataTPED(tpedFilename, hapData->nloci, hapData->nhaps, USE_PMAP);
        }
        else if (VCF) {
            hapData = readHaplotypeDataVCF(vcfFilename,UNPHASED);
            if (CALC_XP || CALC_XPNSL)
            {
                hapData2 = readHaplotypeDataVCF(vcfFilename2,UNPHASED);
                if (hapData->nloci != hapData2->nloci)
                {
                    cerr << "ERROR: Haplotypes from " << vcfFilename << " and " << vcfFilename2 << " do not have the same number of loci.\n";
                    return 1;
                }
            }
            if(!CALC_NSL && !CALC_XPNSL && !USE_PMAP) {
                mapData = readMapData(mapFilename, hapData->nloci, USE_PMAP);
            }
            else{//Load physical positions
                mapData = readMapDataVCF(vcfFilename, hapData->nloci);
            }
        }
        else
        {
            hapData = readHaplotypeData(hapFilename,UNPHASED);
            if (CALC_XP || CALC_XPNSL)
            {
                hapData2 = readHaplotypeData(hapFilename2,UNPHASED);
                if (hapData->nloci != hapData2->nloci)
                {
                    cerr << "ERROR: Haplotypes from " << hapFilename << " and " << hapFilename2 << " do not have the same number of loci.\n";
                    return 1;
                }
            }
            mapData = readMapData(mapFilename, hapData->nloci, USE_PMAP);
        }
    }
    catch (...)
    {
        return 1;
    }

    //Check if map is in order
    for (int i = 1; i < mapData->nloci; i++) {
        if ( mapData->physicalPos[i] < mapData->physicalPos[i - 1] ) {
            cerr << "ERROR: Variant physical position must be monotonically increasing.\n";
            cerr << "\t" << mapData->locusName[i] << " " << mapData->physicalPos[i] << " appears after";
            cerr << "\t" << mapData->locusName[i - 1] << " " << mapData->physicalPos[i - 1] << "\n";
            return 1;
        }
        if ( !CALC_NSL && mapData->geneticPos[i] < mapData->geneticPos[i - 1] ) {
            cerr << "ERROR: Variant genetic position must be monotonically increasing.\n";
            cerr << "\t" << mapData->locusName[i] << " " << mapData->geneticPos[i] << " appears after";
            cerr << "\t" << mapData->locusName[i - 1] << " " << mapData->geneticPos[i - 1] << "\n";
            return 1;
        }
    }


    if (EHH1K >= hapData->nhaps)
    {

    }

    if (SINGLE_EHH)
    {
        queryLoc = queryFound(mapData, query);
        //double queryFreq = calcFreq(hapData, queryLoc);
        if (queryLoc < 0)
        {
            cerr << "ERROR: Could not find specific locus query, " << query << ", in data.\n";
            return 1;
        }
        /*
        else if (queryFreq < MAF || 1 - queryFreq < MAF)
        {
            cerr << "ERROR: EHH for '1' and '0' haplotypes not calculated for " << query << ". MAF < " << MAF << ".\n";
            return 1;
        }
        */
        else
        {
            cerr << "Found " << query << " in data.\n";
        }
    }

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
/*
        MapData *newMapData;
        HaplotypeData *newHapData;
        HaplotypeData *newHapData2;
        double *newfreq1;
        double *newfreq2;
*/
        //int count = 0;
        for (int i = 0; i < hapData->nloci; i++)
        {
            freq1[i] = calcFreq(hapData, i, UNPHASED);
            freq2[i] = calcFreq(hapData2, i, UNPHASED);
            /*
            if(freq1[i] > MAF && 1 - freq1[i] > MAF && 
                freq2[i] > MAF && 1 - freq2[i] > MAF){
                count++;
            }
            */
        }

        //Filtering kills power for XP stats
/*
        if (SKIP) //prefilter all sites < MAF
        {
            cerr << "Removing all variants < " << MAF << " in both pops.\n";
            flog << "Removing all variants < " << MAF << " in both pops.\n";
            newfreq1 = new double [count];
            newfreq2 = new double [count];
            newMapData = initMapData(count);
            newMapData->chr = mapData->chr;
            int j = 0;
            for (int locus = 0; locus < mapData->nloci; locus++)
            {
                if (freq1[locus] > MAF && 1 - freq1[locus] > MAF && 
                    freq2[locus] > MAF && 1 - freq2[locus] > MAF)
                {
                    newMapData->physicalPos[j] = mapData->physicalPos[locus];
                    newMapData->geneticPos[j] = mapData->geneticPos[locus];
                    newMapData->locusName[j] = mapData->locusName[locus];
                    newfreq1[j] = freq1[locus];
                    newfreq2[j] = freq2[locus];
                    j++;
                }
            }

            newHapData = initHaplotypeData(hapData->nhaps, count);
            newHapData2 = initHaplotypeData(hapData2->nhaps, count);

            for (int hap = 0; hap < newHapData->nhaps; hap++)
            {
                j = 0;
                for (int locus = 0; locus < mapData->nloci; locus++)
                {
                    if (freq1[locus] > MAF && 1 - freq1[locus] > MAF && 
                        freq2[locus] > MAF && 1 - freq2[locus] > MAF)
                    {
                        newHapData->data[hap][j] = hapData->data[hap][locus];
                        j++;
                    }
                }
            }

            for (int hap = 0; hap < newHapData2->nhaps; hap++)
            {
                j = 0;
                for (int locus = 0; locus < mapData->nloci; locus++)
                {
                    if (freq1[locus] > MAF && 1 - freq1[locus] > MAF && 
                        freq2[locus] > MAF && 1 - freq2[locus] > MAF)
                    {
                        newHapData2->data[hap][j] = hapData2->data[hap][locus];
                        j++;
                    }
                }
            }

            cerr << "Removed " << mapData->nloci - count << " low frequency variants.\n";
            flog << "Removed " << mapData->nloci - count << " low frequency variants.\n";

            delete [] freq1;
            delete [] freq2;
            freq1 = newfreq1;
            freq2 = newfreq2;
            newfreq1 = NULL;
            newfreq2 = NULL;

            releaseHapData(hapData);
            releaseHapData(hapData2);
            hapData = newHapData;
            hapData2 = newHapData2;
            newHapData = NULL;
            newHapData2 = NULL;

            releaseMapData(mapData);
            mapData = newMapData;
            newMapData = NULL;
        }
*/
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

void calculatePi(HaplotypeData *hapData, MapData *mapData, int winsize, string outFilename)
{
    ofstream fout;
    fout.open(outFilename.c_str());
    if (fout.fail())
    {
        cerr << "ERROR: Failed to open " << outFilename << " for writing.\n";
        throw 1;
    }

    int start = 1;
    //int end = mapData->physicalPos[mapData->nloci - 1];
    int startLocus = 0;
    int endLocus = 0;
    //Identify the start and end indicies in the window
    int winStart = start;
    int winEnd = winStart + winsize - 1;
    int pos;
    int length = 0;
    double pi = 0;
    double denominator = (hapData->nhaps) * (hapData->nhaps - 1) * 0.5;
    int locus;
    for (locus = 0; locus < mapData->nloci; locus++)
    {
        pos = mapData->physicalPos[locus];
        while (pos > winEnd)
        {
            endLocus = locus - 1;
            length = endLocus - startLocus + 1;
            fout << winStart << " " << winEnd << " ";
            //do stuff

            if (length == 0)
            {
                pi = 0;
            }
            else
            {
                for (int i = 0; i < hapData->nhaps; i++)
                {
                    for (int j = i + 1; j < hapData->nhaps; j++)
                    {
                        pi += hamming_dist_ptr(hapData->data[i] + startLocus, hapData->data[j] + startLocus, length);
                    }
                }
            }

            fout << pi / denominator << endl;

            winStart += winsize;
            winEnd += winsize;

            startLocus = locus;
            pi = 0;
        }
    }
    //final window
    endLocus = locus - 1;
    length = endLocus - startLocus + 1;
    fout << winStart << " " << winEnd << " ";
    //do stuff

    if (length == 0)
    {
        pi = 0;
    }
    else
    {
        for (int i = 0; i < hapData->nhaps; i++)
        {
            for (int j = i + 1; j < hapData->nhaps; j++)
            {
                pi += hamming_dist_ptr(hapData->data[i] + startLocus, hapData->data[j] + startLocus, length);
            }
        }
    }

    fout << pi / denominator << endl;



    fout.close();
    return;
}

int queryFound(MapData *mapData, string query)
{
    for (int locus = 0; locus < mapData->nloci; locus++)
    {
        if (mapData->locusName[locus].compare(query) == 0) return locus;
    }

    return -1;
}

double calcFreq(HaplotypeData *hapData, int locus, bool unphased)
{
    double total = 0;
    double freq = 0;

    for (int hap = 0; hap < hapData->nhaps; hap++)
    {
        if (hapData->data[hap][locus] != MISSING_CHAR)
        {
            if (hapData->data[hap][locus] == '1') freq += 1;
            else if (hapData->data[hap][locus] == '2') freq += 2;
            
            if (unphased) total+=2;
            else total++;
        }
    }
    return (freq / total);
}

void query_locus(void *order)
{
    work_order_t *p = (work_order_t *)order;
    char **data = p->hapData->data;
    int nloci = p->hapData->nloci;
    int nhaps = p->hapData->nhaps;
    int *physicalPos = p->mapData->physicalPos;
    double *geneticPos = p->mapData->geneticPos;
    //string *locusName = p->mapData->locusName;
    //ofstream *flog = p->flog;
    ofstream *fout = p->fout;
    string outFilename = p->filename;
    //int SCALE_PARAMETER = p->params->getIntFlag(ARG_GAP_SCALE);
    //int MAX_GAP = p->params->getIntFlag(ARG_MAX_GAP);
    //double EHH_CUTOFF = p->params->getDoubleFlag(ARG_CUTOFF);
    bool ALT = p->params->getBoolFlag(ARG_ALT);
    bool WAGH = p->params->getBoolFlag(ARG_WAGH);
    double (*calc)(map<string, int> &, int, bool) = p->calc;

    int locus = p->queryLoc;
    int queryPad = p->params->getIntFlag(ARG_QWIN);
    int stopLeft = locus;
    for (int i = locus - 1; i >= 0; i--)
    {
        if (physicalPos[locus] - physicalPos[i] <= queryPad) stopLeft = i;
    }
    int stopRight = locus;
    for (int i = locus + 1; i < nloci; i++)
    {
        if (physicalPos[i] - physicalPos[locus] <= queryPad) stopRight = i;
    }

    //EHH to the left of the core snp
    double current_derived_ehh = 1;
    double current_ancestral_ehh = 1;
    double current_ehh = 1;
    double derivedCount = 0;

    //A list of all the haplotypes
    //Starts with just the focal snp and grows outward
    string *haplotypeList = new string[nhaps];
    for (int hap = 0; hap < nhaps; hap++)
    {
        derivedCount += ( data[hap][locus] == '1' ) ? 1 : 0;
        //char digit[2];
        //sprintf(digit, "%d", data[hap][locus]);
        haplotypeList[hap] = data[hap][locus];
    }

    current_ehh = (derivedCount > 1) ? nCk(derivedCount, 2) / nCk(nhaps, 2) : 0;

    /*
        if (derivedCount == 0 || derivedCount == nhaps)
        {
            cerr << "ERROR: " << locusName[locus] << " is monomorphic.\n";
            (*fout) << "ERROR: " << locusName[locus] << " is monomorphic.\n";
            return;
        }
    */
    //cerr << "numHaps: " << nhaps << "\nderivedCounts: " << derivedCount << endl;

    int **ancestralHapColor = new int *[int(nhaps - derivedCount)];
    for (int i = 0; i < nhaps - derivedCount; i++)
    {
        ancestralHapColor[i] = new int[stopRight - stopLeft + 1];
        ancestralHapColor[i][locus - stopLeft] = 0;
    }
    int **derivedHapColor = new int *[int(derivedCount)];
    for (int i = 0; i < derivedCount; i++)
    {
        derivedHapColor[i] = new int[stopRight - stopLeft + 1];
        derivedHapColor[i][locus - stopLeft] = 0;
    }

    //cerr << "allocated hap color arrays.\n";

    string *tempResults = new string[locus - stopLeft];
    int tempIndex = locus - stopLeft - 1;
    int derivedCurrentColor = 0;
    int ancestralCurrentColor = 0;

    for (int i = locus - 1; i >= stopLeft; i--)
    {
        int numDerived = 0;
        int numAncestral = 0;
        map<string, int> ancestralHapCount;
        map<string, int> derivedHapCount;
        map<string, int> hapCount;

        for (int hap = 0; hap < nhaps; hap++)
        {
            bool isDerived = ( data[hap][locus] == '1') ? 1 : 0;
            //build haplotype string
            //char digit[2];
            //sprintf(digit, "%d", data[hap][i]);
            haplotypeList[hap] += data[hap][i];
            string hapStr = haplotypeList[hap];

            if (isDerived)
            {
                //count derived hapoltype
                if (derivedHapCount.count(hapStr) == 0) derivedHapCount[hapStr] = 1;
                else derivedHapCount[hapStr]++;
                numDerived++;
            }
            else
            {
                //count ancestral haplotype
                if (ancestralHapCount.count(hapStr) == 0) ancestralHapCount[hapStr] = 1;
                else ancestralHapCount[hapStr]++;
                numAncestral++;
            }

            if (hapCount.count(hapStr) == 0) hapCount[hapStr] = 1;
            else hapCount[hapStr]++;

        }
        //Write functions to fill in haplotype colors here
        fillColors(derivedHapColor, derivedHapCount, haplotypeList, nhaps, tempIndex, derivedCurrentColor, true);
        fillColors(ancestralHapColor, ancestralHapCount, haplotypeList, nhaps, tempIndex, ancestralCurrentColor, true);

        current_derived_ehh = (*calc)(derivedHapCount, numDerived, ALT);
        current_ancestral_ehh = (*calc)(ancestralHapCount, numAncestral, ALT);
        current_ehh = (*calc)(hapCount, numAncestral + numDerived, ALT);
        char tempStr[100];
        sprintf(tempStr, "%d\t%f\t%f\t%f\t%f", physicalPos[i] - physicalPos[locus], geneticPos[i] - geneticPos[locus], current_derived_ehh, current_ancestral_ehh, current_ehh);
        tempResults[tempIndex] = string(tempStr);
        tempIndex--;
    }

    delete [] haplotypeList;

    for (int i = 0; i < locus - stopLeft; i++)
    {
        (*fout) << tempResults[i] << "\n";
    }
    delete [] tempResults;

    //calculate EHH to the right
    current_derived_ehh = 1;
    current_ancestral_ehh = 1;
    current_ehh = (derivedCount > 1) ? nCk(derivedCount, 2) / nCk(nhaps, 2) : 0;

    fout->precision(6);
    (*fout) << std::fixed <<  physicalPos[locus] - physicalPos[locus]  << "\t"
            << geneticPos[locus] - geneticPos[locus] << "\t"
            << current_derived_ehh << "\t"
            << current_ancestral_ehh << "\t"
            << current_ehh << endl;

    //A list of all the haplotypes
    //Starts with just the focal snp and grows outward
    haplotypeList = new string[nhaps];
    for (int hap = 0; hap < nhaps; hap++)
    {
        //char digit[2];
        //sprintf(digit, "%d", data[hap][locus]);
        haplotypeList[hap] = data[hap][locus];
    }

    derivedCurrentColor = 0;
    ancestralCurrentColor = 0;

    //while(current_ancestral_ehh > EHH_CUTOFF || current_derived_ehh > EHH_CUTOFF)
    for (int i = locus + 1; i <= stopRight; i++)
    {
        int numDerived = 0;
        int numAncestral = 0;
        map<string, int> ancestralHapCount;
        map<string, int> derivedHapCount;
        map<string, int> hapCount;

        for (int hap = 0; hap < nhaps; hap++)
        {
            bool isDerived = ( data[hap][locus] == '1' ) ? 1 : 0;
            //build haplotype string
            //char digit[2];
            //sprintf(digit, "%d", data[hap][i]);
            haplotypeList[hap] += data[hap][i];
            string hapStr = haplotypeList[hap];

            if (isDerived)
            {
                //count hapoltypes
                if (derivedHapCount.count(hapStr) == 0) derivedHapCount[hapStr] = 1;
                else derivedHapCount[hapStr]++;
                numDerived++;
            }
            else //ancestral
            {
                if (ancestralHapCount.count(hapStr) == 0) ancestralHapCount[hapStr] = 1;
                else ancestralHapCount[hapStr]++;
                numAncestral++;
            }

            if (hapCount.count(hapStr) == 0) hapCount[hapStr] = 1;
            else hapCount[hapStr]++;
        }

        //Write functions to fill in haplotype colors here
        fillColors(derivedHapColor, derivedHapCount, haplotypeList, nhaps, i - stopLeft, derivedCurrentColor, false);
        fillColors(ancestralHapColor, ancestralHapCount, haplotypeList, nhaps, i - stopLeft, ancestralCurrentColor, false);

        current_derived_ehh = (*calc)(derivedHapCount, numDerived, ALT);
        current_ancestral_ehh = (*calc)(ancestralHapCount, numAncestral, ALT);
        current_ehh = (*calc)(hapCount, numAncestral + numDerived, ALT);

        (*fout) << physicalPos[i] - physicalPos[locus] << "\t"
                << geneticPos[i] - geneticPos[locus] << "\t"
                << current_derived_ehh << "\t"
                << current_ancestral_ehh << "\t"
                << current_ehh << endl;
    }

    delete [] haplotypeList;

    ofstream fout2;
    string outFilenameDer = outFilename + ".der.colormap";
    fout2.open(outFilenameDer.c_str());
    if (fout2.fail())
    {
        cerr << "ERROR: could not open " << outFilenameDer << " for writing.\n";
        throw 1;
    }

    for (int i = 0; i < derivedCount; i++)
    {
        for (int j = 0; j < stopRight - stopLeft + 1; j++)
        {
            fout2 << derivedHapColor[i][j] << " ";
        }
        fout2 << endl;
    }
    fout2.close();

    string outFilenameAnc = outFilename + ".anc.colormap";
    fout2.open(outFilenameAnc.c_str());
    if (fout2.fail())
    {
        cerr << "ERROR: could not open " << outFilenameAnc << " for writing.\n";
        throw 1;
    }

    for (int i = 0; i < nhaps - derivedCount; i++)
    {
        for (int j = 0; j < stopRight - stopLeft + 1; j++)
        {
            fout2 << ancestralHapColor[i][j] << " ";
        }
        fout2 << endl;
    }
    fout2.close();
    return;
}

void query_locus_soft(void *order)
{
    work_order_t *p = (work_order_t *)order;
    char **data = p->hapData->data;
    int nloci = p->hapData->nloci;
    int nhaps = p->hapData->nhaps;
    int *physicalPos = p->mapData->physicalPos;
    double *geneticPos = p->mapData->geneticPos;
    //string *locusName = p->mapData->locusName;

    //ofstream *flog = p->flog;
    ofstream *fout = p->fout;
    string outFilename = p->filename;
    //int SCALE_PARAMETER = p->params->getIntFlag(ARG_GAP_SCALE);
    //int MAX_GAP = p->params->getIntFlag(ARG_MAX_GAP);
    //double EHH_CUTOFF = p->params->getDoubleFlag(ARG_CUTOFF);
    //bool ALT = p->params->getBoolFlag(ARG_ALT);
    //bool WAGH = p->params->getBoolFlag(ARG_WAGH);

    int locus = p->queryLoc;
    int queryPad = p->params->getIntFlag(ARG_QWIN);
    int stopLeft = locus;
    for (int i = locus - 1; i >= 0; i--)
    {
        if (physicalPos[locus] - physicalPos[i] <= queryPad) stopLeft = i;
    }
    int stopRight = locus;
    for (int i = locus + 1; i < nloci; i++)
    {
        if (physicalPos[i] - physicalPos[locus] <= queryPad) stopRight = i;
    }

    //EHH to the left of the core snp

    double current_ehh1 = 1;
    double current_ehh2d1 = 1;
    double current_ehh12 = 1;

    double previous_ehh1 = 1;
    double previous_ehh2d1 = 1;
    double previous_ehh12 = 1;

    double derivedCount = 0;
    //A list of all the haplotypes
    //Starts with just the focal snp and grows outward
    map<string, int> tempHapCount;
    string *haplotypeList = new string[nhaps];
    for (int hap = 0; hap < nhaps; hap++)
    {
        derivedCount += data[hap][locus];
        char digit[2];
        sprintf(digit, "%d", data[hap][locus]);
        haplotypeList[hap] = digit;
        string hapStr = haplotypeList[hap];
        //count hapoltype freqs
        if (tempHapCount.count(hapStr) == 0) tempHapCount[hapStr] = 1;
        else tempHapCount[hapStr]++;
    }

    triplet_t res = calculateSoft(tempHapCount, nhaps);
    current_ehh1 = res.h1;
    current_ehh2d1 = res.h2dh1;
    current_ehh12 = res.h12;
    previous_ehh1 = res.h1;
    previous_ehh2d1 = res.h2dh1;
    previous_ehh12 = res.h12;

    //cerr << "numHaps: " << nhaps << "\nderivedCounts: " << derivedCount << endl;
    /*
    int **ancestralHapColor = new int*[int(nhaps-derivedCount)];
    for(int i = 0; i < nhaps-derivedCount; i++)
    {
      ancestralHapColor[i] = new int[stopRight-stopLeft+1];
      ancestralHapColor[i][locus-stopLeft] = 0;
    }
    int **derivedHapColor = new int*[int(derivedCount)];
    for(int i = 0; i < derivedCount; i++)
    {
      derivedHapColor[i] = new int[stopRight-stopLeft+1];
      derivedHapColor[i][locus-stopLeft] = 0;
    }
    */
    //cerr << "allocated hap color arrays.\n";

    string *tempResults = new string[locus - stopLeft];
    int tempIndex = locus - stopLeft - 1;
    //int derivedCurrentColor = 0;
    //int ancestralCurrentColor = 0;

    for (int i = locus - 1; i >= stopLeft; i--)
    {
        //int numDerived = 0;
        //int numAncestral = 0;
        map<string, int> hapCount;

        for (int hap = 0; hap < nhaps; hap++)
        {
            //build haplotype string
            char digit[2];
            sprintf(digit, "%d", data[hap][i]);
            haplotypeList[hap] += digit;
            string hapStr = haplotypeList[hap];

            //count hapoltype freqs
            if (hapCount.count(hapStr) == 0) hapCount[hapStr] = 1;
            else hapCount[hapStr]++;
        }

        //We've now counted all of the unique haplotypes extending out of the core SNP
        res = calculateSoft(hapCount, nhaps);
        current_ehh1 = res.h1;
        current_ehh2d1 = res.h2dh1;
        current_ehh12 = res.h12;
        //Write functions to fill in haplotype colors here
        /*
        fillColors(derivedHapColor, derivedHapCount,haplotypeList, nhaps,tempIndex,derivedCurrentColor,true);
        fillColors(ancestralHapColor, ancestralHapCount,haplotypeList, nhaps,tempIndex,ancestralCurrentColor,true);
        */

        char tempStr[100];
        sprintf(tempStr, "%d\t%f\t%f\t%f\t%f", physicalPos[i] - physicalPos[locus], geneticPos[i] - geneticPos[locus], current_ehh1, current_ehh12, current_ehh2d1);
        tempResults[tempIndex] = string(tempStr);
        tempIndex--;
    }

    delete [] haplotypeList;

    for (int i = 0; i < locus - stopLeft; i++)
    {
        (*fout) << tempResults[i] << "\n";
    }
    delete [] tempResults;

    //calculate EHH to the right
    current_ehh1 = 1;
    current_ehh2d1 = 1;
    current_ehh12 = 1;

    previous_ehh1 = 1;
    previous_ehh2d1 = 1;
    previous_ehh12 = 1;

    //A list of all the haplotypes
    //Starts with just the focal snp and grows outward
    tempHapCount.clear();
    haplotypeList = new string[nhaps];
    for (int hap = 0; hap < nhaps; hap++)
    {
        char digit[2];
        sprintf(digit, "%d", data[hap][locus]);
        haplotypeList[hap] = digit;
        string hapStr = haplotypeList[hap];
        //count hapoltype freqs
        if (tempHapCount.count(hapStr) == 0) tempHapCount[hapStr] = 1;
        else tempHapCount[hapStr]++;
    }

    res = calculateSoft(tempHapCount, nhaps);
    current_ehh1 = res.h1;
    current_ehh2d1 = res.h2dh1;
    current_ehh12 = res.h12;
    previous_ehh1 = res.h1;
    previous_ehh2d1 = res.h2dh1;
    previous_ehh12 = res.h12;

    fout->precision(6);
    (*fout) << std::fixed <<  physicalPos[locus] - physicalPos[locus]  << "\t"
            << geneticPos[locus] - geneticPos[locus] << "\t"
            << current_ehh1 << "\t"
            << current_ehh12 << "\t"
            << current_ehh2d1 << endl;

    //while(current_ancestral_ehh > EHH_CUTOFF || current_derived_ehh > EHH_CUTOFF)
    for (int i = locus + 1; i <= stopRight; i++)
    {
        map<string, int> hapCount;

        for (int hap = 0; hap < nhaps; hap++)
        {
            //build haplotype string
            char digit[2];
            sprintf(digit, "%d", data[hap][i]);
            haplotypeList[hap] += digit;
            string hapStr = haplotypeList[hap];

            //count hapoltype freqs
            if (hapCount.count(hapStr) == 0) hapCount[hapStr] = 1;
            else hapCount[hapStr]++;
        }

        //We've now counted all of the unique haplotypes extending out of the core SNP
        res = calculateSoft(hapCount, nhaps);
        current_ehh1 = res.h1;
        current_ehh2d1 = res.h2dh1;
        current_ehh12 = res.h12;

        //Write functions to fill in haplotype colors here
        /*
        fillColors(derivedHapColor, derivedHapCount,haplotypeList, nhaps,i-stopLeft,derivedCurrentColor,false);
        fillColors(ancestralHapColor, ancestralHapCount,haplotypeList, nhaps,i-stopLeft,ancestralCurrentColor,false);
        */

        (*fout) << physicalPos[i] - physicalPos[locus] << "\t"
                << geneticPos[i] - geneticPos[locus] << "\t"
                << current_ehh1 << "\t"
                << current_ehh12 << "\t"
                << current_ehh2d1 << endl;
    }

    delete [] haplotypeList;

    /*
    ofstream fout2;
    string outFilenameDer = outFilename + ".der.colormap";
    fout2.open(outFilenameDer.c_str());
    if(fout2.fail())
    {
      cerr << "ERROR: could not open " << outFilenameDer << " for writing.\n";
      throw 1;
    }
    for(int i = 0; i < derivedCount; i++)
    {
      for(int j = 0; j < stopRight-stopLeft+1; j++)
      {
        fout2 << derivedHapColor[i][j] << " ";
        }
      fout2 << endl;
    }
    fout2.close();
    string outFilenameAnc = outFilename + ".anc.colormap";
    fout2.open(outFilenameAnc.c_str());
    if(fout2.fail())
    {
      cerr << "ERROR: could not open " << outFilenameAnc << " for writing.\n";
      throw 1;
    }
    for(int i = 0; i < nhaps-derivedCount; i++)
    {
      for(int j = 0; j < stopRight-stopLeft+1; j++)
      {
        fout2 << ancestralHapColor[i][j] << " ";
        }
      fout2 << endl;
    }
    fout2.close();
    */
    return;
}

void fillColors(int **hapColor, map<string, int> &hapCount, string *haplotypeList, int hapListLength, int currentLoc, int &currentColor, bool left)
{
    map<string, int>::iterator i;
    //int numUniqueHaps = hapCount.size();
    //int mostCommonHapCount = 0;
    int nhaps = 0;

    for (i = hapCount.begin(); i != hapCount.end(); i++)
    {
        nhaps += i->second;
    }
    //holds colors for haplotypes that have already been seen in the search
    map<string, int> hapSeen;

    string mostCommonHap = "_NONE_";

    int colorIndex = 0;
    int previousLoc = (left) ? currentLoc + 1 : currentLoc - 1;
    for (int j = 0; j < hapListLength; j++)
    {
        string hapStr = haplotypeList[j];
        if (hapCount.count(hapStr) == 0) continue;

        /*
        for(int h = 0; h < hapListLength; h++)
        {
          if(h == j) cerr << ">>";
          else cerr << "  ";
          cerr << haplotypeList[h] << endl;
        }
        */

        if (hapCount[hapStr] == 1)
        {
            //cerr << "Unique hap\n";
            hapColor[colorIndex][currentLoc] = -1;
            colorIndex++;
            continue;
        }

        //If there was a split in the current haplotype family
        //let the most common continued haplotype keep the color
        //then increment the color for the less common one
        if (familyDidSplit(hapStr, hapCount[hapStr], hapColor, nhaps, colorIndex, previousLoc, mostCommonHap))
        {
            //cerr << "split\n";
            //cerr << hapStr << " " << mostCommonHap << endl;
            if (hapStr.compare(mostCommonHap) == 0)
            {
                //cerr << "common\n";
                hapColor[colorIndex][currentLoc] = hapColor[colorIndex][previousLoc];
            }
            else
            {
                //cerr << "not common ";
                if (hapSeen.count(hapStr) == 0) //not seen
                {
                    //cerr << "not seen\n";
                    hapColor[colorIndex][currentLoc] = ++currentColor;
                    hapSeen[hapStr] = currentColor;
                }
                else // seen
                {
                    //cerr << "seen\n";
                    hapColor[colorIndex][currentLoc] = hapSeen[hapStr];
                }
            }
        }
        //Family did not split, so keep the color it had
        else
        {
            //cerr << "no split\n";
            hapColor[colorIndex][currentLoc] = hapColor[colorIndex][previousLoc];
        }
        colorIndex++;

        //string junk;
        //cin >> junk;
    }
    return;
}

bool familyDidSplit(const string &hapStr, const int hapCount, int **hapColor, const int nhaps, const int colorIndex, const int previousLoc, string &mostCommonHap)
{
    //cerr << "most common hap: " << mostCommonHap << endl;
    int previousColor = hapColor[colorIndex][previousLoc];
    int numPrevColor = 0;

    for (int i = 0; i < nhaps; i++)
    {
        if (hapColor[i][previousLoc] == previousColor) numPrevColor++;
    }

    //cerr << "numPrevColor: " << numPrevColor << "\nhapCount: " << hapCount << endl;

    if (numPrevColor == hapCount) return false;

    if ( hapCount > double(numPrevColor) / 2.0  )
    {
        mostCommonHap = hapStr;
        return true;
    }
    else if (mostCommonHap.compare("_NONE_") == 0 && hapCount == double(numPrevColor) / 2.0)
    {
        mostCommonHap = hapStr;
        return true;
    }
    else return true;
}

void calc_ihs(void *order)
{
    work_order_t *p = (work_order_t *)order;
    char **data = p->hapData->data;
    int nloci = p->hapData->nloci;
    int nhaps = p->hapData->nhaps;
    int *physicalPos = p->mapData->physicalPos;
    double *geneticPos = p->mapData->geneticPos;
    string *locusName = p->mapData->locusName;
    int id = p->id;
    double *ihs = p->ihs;
    double *ihh1 = p->ihh1;
    double *ihh2 = p->ihh2;
    double *ihhDerivedLeft    = p->ihhDerivedLeft;
    double *ihhDerivedRight   = p->ihhDerivedRight;
    double *ihhAncestralLeft  = p->ihhAncestralLeft;
    double *ihhAncestralRight = p->ihhAncestralRight;
    double *freq = p->freq;
    ofstream *flog = p->flog;
    Bar *pbar = p->bar;

    int SCALE_PARAMETER = p->params->getIntFlag(ARG_GAP_SCALE);
    int MAX_GAP = p->params->getIntFlag(ARG_MAX_GAP);
    double EHH_CUTOFF = p->params->getDoubleFlag(ARG_CUTOFF);
    bool ALT = p->params->getBoolFlag(ARG_ALT);
    bool WAGH = p->params->getBoolFlag(ARG_WAGH);
    bool TRUNC = p->params->getBoolFlag(ARG_TRUNC);
    double MAF = p->params->getDoubleFlag(ARG_MAF);
    int numThreads = p->params->getIntFlag(ARG_THREAD);
    int MAX_EXTEND = ( p->params->getIntFlag(ARG_MAX_EXTEND) <= 0 ) ? physicalPos[nloci - 1] - physicalPos[0] : p->params->getIntFlag(ARG_MAX_EXTEND);
    //bool SKIP = p->params->getBoolFlag(ARG_SKIP);
    bool WRITE_DETAILED_IHS = p->params->getBoolFlag(ARG_IHS_DETAILED);
    double (*calc)(map<string, int> &, int, bool) = p->calc;

    bool unphased = p->params->getBoolFlag(ARG_UNPHASED);

    //int step = (stop - start) / (pbar->totalTicks);
    int step = (nloci / numThreads) / (pbar->totalTicks);
    if (step == 0) step = 1;

    bool isDerived;
    string hapStr;

    for (int locus = id; locus < nloci; locus += numThreads)
    {
        if (locus % step == 0) advanceBar(*pbar, double(step));

        ihs[locus] = 0;
        //freq[locus] = MISSING;
        ihh1[locus] = MISSING;
        ihh2[locus] = MISSING;
        if (WRITE_DETAILED_IHS) {
            ihhDerivedLeft[locus]    = MISSING;
            ihhDerivedRight[locus]   = MISSING;
            ihhAncestralLeft[locus]  = MISSING;
            ihhAncestralRight[locus] = MISSING;
        }
        bool skipLocus = 0;
        //If the focal snp has MAF < MAF, then skip this locus
        if (freq[locus] < MAF || freq[locus] > 1 - MAF)
        {
            pthread_mutex_lock(&mutex_log);
            (*flog) << "WARNING: Locus " << locusName[locus]
                    << " has MAF < " << MAF << ". Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus] << "\n";
            pthread_mutex_unlock(&mutex_log);
            ihs[locus] = MISSING;
            skipLocus = 0;
            continue;
        }

        //EHH to the left of the core snp
        double current_derived_ehh = 1;
        double current_ancestral_ehh = 1;
        double previous_derived_ehh = 1;
        double previous_ancestral_ehh = 1;
        int currentLocus = locus;
        int nextLocus = locus - 1;
        double derived_ihh = 0;
        double ancestral_ihh = 0;

        double derived_ihh_left    = 0;
        double ancestral_ihh_left  = 0;
        double derived_ihh_right   = 0;
        double ancestral_ihh_right = 0;

        //used for unphased analyses
        double current_notDerived_ehh = 1;
        double current_notAncestral_ehh = 1;
        double previous_notDerived_ehh = 1;
        double previous_notAncestral_ehh = 1;
        double notDerived_ihh = 0;
        double notAncestral_ihh = 0;

        //left right for unphased
        double notDerived_ihh_left    = 0;
        double notAncestral_ihh_left  = 0;
        double notDerived_ihh_right   = 0;
        double notAncestral_ihh_right = 0;

        //A list of all the haplotypes
        //Starts with just the focal snp and grows outward
        string *haplotypeList = new string[nhaps];
        for (int hap = 0; hap < nhaps; hap++)
        {
            //derivedCount += ( data[hap][locus] == '1' ) ? 1 : 0;
            haplotypeList[hap] = data[hap][locus];
        }

        while (current_derived_ehh > EHH_CUTOFF || current_ancestral_ehh > EHH_CUTOFF)
        {
            if (nextLocus < 0)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << EHH_CUTOFF
                        << ". ";
                if (!TRUNC){
                    skipLocus = 1;
                    (*flog) << "Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus];
                }
                (*flog) << "\n";
                pthread_mutex_unlock(&mutex_log);
                break;
            }
            else if (physicalPos[currentLocus] - physicalPos[nextLocus] > MAX_GAP)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached a gap of " << physicalPos[currentLocus] - physicalPos[nextLocus]
                        << "bp > " << MAX_GAP << "bp. Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus] << "\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }

            //Check to see if the gap between the markers is huge, if so, scale it in an ad hoc way as in
            //Voight et al. (2006)
            double scale = double(SCALE_PARAMETER) / double(physicalPos[currentLocus] - physicalPos[nextLocus]);
            if (scale > 1) scale = 1;

            currentLocus = nextLocus;
            nextLocus--;

            int numDerived = 0;
            int numAncestral = 0;
            map<string, int> ancestralHapCount;
            map<string, int> derivedHapCount;

            //used for unphased analyses
            int numHet = 0;
            map<string, int> notAncestralHapCount;
            map<string, int> notDerivedHapCount;

            for (int hap = 0; hap < nhaps; hap++)
            {
                if(unphased){
                    haplotypeList[hap] += data[hap][currentLocus];
                    hapStr = haplotypeList[hap];
                    if (data[hap][locus] == '0'){
                        //count ancestral haplotype
                        if (ancestralHapCount.count(hapStr) == 0) ancestralHapCount[hapStr] = 1;
                        else ancestralHapCount[hapStr]++;
                        numAncestral++;
                        //count non-derived hapoltype
                        if (notDerivedHapCount.count(hapStr) == 0) notDerivedHapCount[hapStr] = 1;
                        else notDerivedHapCount[hapStr]++;
                    }
                    else if (data[hap][locus] == '1'){
                        //count non-derived hapoltype
                        if (notDerivedHapCount.count(hapStr) == 0) notDerivedHapCount[hapStr] = 1;
                        else notDerivedHapCount[hapStr]++;
                        //count non-ancestral haplotype
                        if (notAncestralHapCount.count(hapStr) == 0) notAncestralHapCount[hapStr] = 1;
                        else notAncestralHapCount[hapStr]++;
                        numHet++;
                    }
                    else{//data[hap][locus] == '2'
                        //count derived hapoltype
                        if (derivedHapCount.count(hapStr) == 0) derivedHapCount[hapStr] = 1;
                        else derivedHapCount[hapStr]++;
                        numDerived++;
                        //count non-ancestral haplotype
                        if (notAncestralHapCount.count(hapStr) == 0) notAncestralHapCount[hapStr] = 1;
                        else notAncestralHapCount[hapStr]++;
                    }
                }
                else{
                    isDerived = ( data[hap][locus] == '1' ) ? 1 : 0;
                    haplotypeList[hap] += data[hap][currentLocus];
                    hapStr = haplotypeList[hap];

                    if (isDerived)
                    {
                        //count derived hapoltype
                        if (derivedHapCount.count(hapStr) == 0) derivedHapCount[hapStr] = 1;
                        else derivedHapCount[hapStr]++;
                        numDerived++;
                    }
                    else
                    {
                        //count ancestral haplotype
                        if (ancestralHapCount.count(hapStr) == 0) ancestralHapCount[hapStr] = 1;
                        else ancestralHapCount[hapStr]++;
                        numAncestral++;
                    }
                }
            }

            //We've now counted all of the unique haplotypes extending out of the core SNP
            //If locus is monomorphic, shoot a warning and skip locus
            //This probably isn't necessary any more
            if ( !unphased && (numDerived == 0 || numAncestral == 0) ) 
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: locus " << locusName[locus]
                        << " (number " << locus + 1 << ") is monomorphic. Skipping calculation at this locus.\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }

            double freqHetGT = double(numHet) / double(numDerived + numAncestral + numHet);
            //double freqAncestralGT = double(numAncestral) / double(numDerived + numAncestral + numHet);

            if ( unphased && freqHetGT > 1-MAF ) 
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: locus " << locusName[locus]
                        << " (number " << locus + 1 << ") has too many hets. Skipping calculation at this locus. "
                        << "het: " << numHet << " hom0: " << numAncestral << " hom1: " << numDerived << ".\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }

            if (current_derived_ehh > EHH_CUTOFF)
            {
                current_derived_ehh = (*calc)(derivedHapCount, numDerived, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                derived_ihh += 0.5 * scale * (geneticPos[currentLocus + 1] - geneticPos[currentLocus]) * (current_derived_ehh + previous_derived_ehh);
                previous_derived_ehh = current_derived_ehh;
            }

            if (current_ancestral_ehh > EHH_CUTOFF)
            {
                current_ancestral_ehh = (*calc)(ancestralHapCount, numAncestral, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                ancestral_ihh += 0.5 * scale * (geneticPos[currentLocus + 1] - geneticPos[currentLocus]) * (current_ancestral_ehh + previous_ancestral_ehh);
                previous_ancestral_ehh = current_ancestral_ehh;
            }

            if(unphased){
                current_notDerived_ehh = (*calc)(notDerivedHapCount, numAncestral + numHet, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                notDerived_ihh += 0.5 * scale * (current_notDerived_ehh + previous_notDerived_ehh);
                previous_notDerived_ehh = current_notDerived_ehh;

                current_notAncestral_ehh = (*calc)(notAncestralHapCount, numDerived + numHet, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                notAncestral_ihh += 0.5 * scale * (current_notAncestral_ehh + previous_notAncestral_ehh);
                previous_notAncestral_ehh = current_notAncestral_ehh;
            }

            //check if currentLocus is beyond MAX_EXTEND
            if (physicalPos[locus] - physicalPos[currentLocus] >= MAX_EXTEND) break;
        }

        delete [] haplotypeList;

        if (skipLocus == 1)
        {
            ihs[locus] = MISSING;
            skipLocus = 0;
            continue;
        }

        derived_ihh_left   = derived_ihh;
        ancestral_ihh_left = ancestral_ihh;

        if(unphased){
            notDerived_ihh_left    = notDerived_ihh;    
            notAncestral_ihh_left  = notAncestral_ihh;
        }

        if(unphased){
                ihh1[locus] = log10(derived_ihh / notDerived_ihh);
                ihh2[locus] = log10(ancestral_ihh / notAncestral_ihh);
                ihs[locus] = (ihh1[locus] > ihh2[locus]) ? ihh1[locus] : 0-ihh2[locus];
            }

        //calculate EHH to the right
        current_derived_ehh = 1;
        current_ancestral_ehh = 1;
        previous_derived_ehh = 1;
        previous_ancestral_ehh = 1;
        currentLocus = locus;
        nextLocus = locus + 1;
        skipLocus = 0;

        //used for unphased analyses
        current_notDerived_ehh = 1;
        current_notAncestral_ehh = 1;
        previous_notDerived_ehh = 1;
        previous_notAncestral_ehh = 1;
        notDerived_ihh = 0;
        notAncestral_ihh = 0;

        //A list of all the haplotypes
        //Starts with just the focal snp and grows outward
        haplotypeList = new string[nhaps];
        for (int hap = 0; hap < nhaps; hap++)
        {
            //char digit[2];
            //sprintf(digit, "%d", data[hap][locus]);
            haplotypeList[hap] = data[hap][locus];
        }

        while (current_ancestral_ehh > EHH_CUTOFF || current_derived_ehh > EHH_CUTOFF)
        {
            if (nextLocus > nloci - 1)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << EHH_CUTOFF
                        << ". ";
                if (!TRUNC){
                    skipLocus = 1;
                    (*flog) << "Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus];
                }
                (*flog) << "\n";
                pthread_mutex_unlock(&mutex_log);
                break;
            }
            else if (physicalPos[nextLocus] - physicalPos[currentLocus] > MAX_GAP)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached a gap of " << physicalPos[nextLocus] - physicalPos[currentLocus]
                        << "bp > " << MAX_GAP << "bp. Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus] << "\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }

            double scale = double(SCALE_PARAMETER) / double(physicalPos[nextLocus] - physicalPos[currentLocus]);
            if (scale > 1) scale = 1;

            currentLocus = nextLocus;
            nextLocus++;

            int numDerived = 0;
            int numAncestral = 0;
            map<string, int> ancestralHapCount;
            map<string, int> derivedHapCount;

            //used for unphased analyses
            int numHet = 0;
            map<string, int> notAncestralHapCount;
            map<string, int> notDerivedHapCount;

            for (int hap = 0; hap < nhaps; hap++)
            {
                if(unphased){
                    haplotypeList[hap] += data[hap][currentLocus];
                    hapStr = haplotypeList[hap];
                    if (data[hap][locus] == '0'){
                        //count ancestral haplotype
                        if (ancestralHapCount.count(hapStr) == 0) ancestralHapCount[hapStr] = 1;
                        else ancestralHapCount[hapStr]++;
                        numAncestral++;
                        //count non-derived hapoltype
                        if (notDerivedHapCount.count(hapStr) == 0) notDerivedHapCount[hapStr] = 1;
                        else notDerivedHapCount[hapStr]++;
                    }
                    else if (data[hap][locus] == '1'){
                        //count non-derived hapoltype
                        if (notDerivedHapCount.count(hapStr) == 0) notDerivedHapCount[hapStr] = 1;
                        else notDerivedHapCount[hapStr]++;
                        //count non-ancestral haplotype
                        if (notAncestralHapCount.count(hapStr) == 0) notAncestralHapCount[hapStr] = 1;
                        else notAncestralHapCount[hapStr]++;
                        numHet++;
                    }
                    else{//data[hap][locus] == '2'
                        //count derived hapoltype
                        if (derivedHapCount.count(hapStr) == 0) derivedHapCount[hapStr] = 1;
                        else derivedHapCount[hapStr]++;
                        numDerived++;
                        //count non-ancestral haplotype
                        if (notAncestralHapCount.count(hapStr) == 0) notAncestralHapCount[hapStr] = 1;
                        else notAncestralHapCount[hapStr]++;
                    }
                }
                else{
                    isDerived = ( data[hap][locus] == '1') ? 1 : 0;
                    haplotypeList[hap] += data[hap][currentLocus];
                    hapStr = haplotypeList[hap];

                    if (isDerived)
                    {
                        //count hapoltypes
                        if (derivedHapCount.count(hapStr) == 0) derivedHapCount[hapStr] = 1;
                        else derivedHapCount[hapStr]++;
                        numDerived++;
                    }
                    else //ancestral
                    {
                        if (ancestralHapCount.count(hapStr) == 0) ancestralHapCount[hapStr] = 1;
                        else ancestralHapCount[hapStr]++;
                        numAncestral++;
                    }
                }
            }

            //We've now counted all of the unique haplotypes extending out of the core SNP
            //If there are no derived alleles at a locus, shoot a warning and skip locus
            if (numDerived == 0 || numAncestral == 0)
            {
                //(*flog) << "WARNING: locus " << locusName[locus]
                //   << " (number " << locus+1 << ") is monomorphic.  Skipping calculation at this locus.\n";
                skipLocus = 1;
                break;
            }

            if (current_derived_ehh > EHH_CUTOFF)
            {
                current_derived_ehh = (*calc)(derivedHapCount, numDerived, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                derived_ihh += 0.5 * scale * (geneticPos[currentLocus] - geneticPos[currentLocus - 1]) * (current_derived_ehh + previous_derived_ehh);
                previous_derived_ehh = current_derived_ehh;
            }

            if (current_ancestral_ehh > EHH_CUTOFF)
            {
                current_ancestral_ehh = (*calc)(ancestralHapCount, numAncestral, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                ancestral_ihh += 0.5 * scale * (geneticPos[currentLocus] - geneticPos[currentLocus - 1]) * (current_ancestral_ehh + previous_ancestral_ehh);
                previous_ancestral_ehh = current_ancestral_ehh;
            }

            if(unphased){
                current_notDerived_ehh = (*calc)(notDerivedHapCount, numAncestral + numHet, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                notDerived_ihh += 0.5 * scale * (current_notDerived_ehh + previous_notDerived_ehh);
                previous_notDerived_ehh = current_notDerived_ehh;

                current_notAncestral_ehh = (*calc)(notAncestralHapCount, numDerived + numHet, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                notAncestral_ihh += 0.5 * scale * (current_notAncestral_ehh + previous_notAncestral_ehh);
                previous_notAncestral_ehh = current_notAncestral_ehh;
            }

            //check if currentLocus is beyond 1Mb
            if (physicalPos[currentLocus] - physicalPos[locus] >= MAX_EXTEND) break;

        }

        derived_ihh_right   = derived_ihh   - derived_ihh_left;
        ancestral_ihh_right = ancestral_ihh - ancestral_ihh_left;

        if(unphased){
            notDerived_ihh_right = notDerived_ihh - notDerived_ihh_left;
            notAncestral_ihh_right = notAncestral_ihh - notAncestral_ihh_left;
        }

        delete [] haplotypeList;

        if (skipLocus == 1)
        {
            ihs[locus] = MISSING;
            skipLocus = 0;
            continue;
        }

        if (ihs[locus] != MISSING)
        {
            if(unphased){
                ihh1[locus] = log10(derived_ihh / notDerived_ihh);
                ihh2[locus] = log10(ancestral_ihh / notAncestral_ihh);
                ihs[locus] = (ihh1[locus] > ihh2[locus]) ? ihh1[locus] : 0-ihh2[locus];
            }
            else{
                ihh1[locus] = derived_ihh;
                ihh2[locus] = ancestral_ihh;
                ihs[locus] = log10(derived_ihh / ancestral_ihh);
            }
            if (WRITE_DETAILED_IHS) {
                if(unphased){//for the time being this is going to report same as phased
                    ihhDerivedLeft[locus]    = derived_ihh_left;
                    ihhDerivedRight[locus]   = derived_ihh_right;
                    ihhAncestralLeft[locus]  = ancestral_ihh_left;
                    ihhAncestralRight[locus] = ancestral_ihh_right;
                }
                else{
                    ihhDerivedLeft[locus]    = derived_ihh_left;
                    ihhDerivedRight[locus]   = derived_ihh_right;
                    ihhAncestralLeft[locus]  = ancestral_ihh_left;
                    ihhAncestralRight[locus] = ancestral_ihh_right;
                }
            }
            //freq[locus] = double(derivedCount) / double(nhaps);
        }
    }
}

void calc_nsl(void *order)
{
    work_order_t *p = (work_order_t *)order;
    char **data = p->hapData->data;
    int nloci = p->hapData->nloci;
    int nhaps = p->hapData->nhaps;
    int *physicalPos = p->mapData->physicalPos;
    //double *geneticPos = p->mapData->geneticPos;
    string *locusName = p->mapData->locusName;
    int id = p->id;
    double *ihs = p->ihs;
    double *ihh1 = p->ihh1;
    double *ihh2 = p->ihh2;
    double *freq = p->freq;
    ofstream *flog = p->flog;
    Bar *pbar = p->bar;

    int SCALE_PARAMETER = p->params->getIntFlag(ARG_GAP_SCALE);
    int MAX_GAP = p->params->getIntFlag(ARG_MAX_GAP);
    //double EHH_CUTOFF = p->params->getDoubleFlag(ARG_CUTOFF);
    double EHH_CUTOFF = 0;
    bool ALT = p->params->getBoolFlag(ARG_ALT);
    bool WAGH = p->params->getBoolFlag(ARG_WAGH);
    bool TRUNC = p->params->getBoolFlag(ARG_TRUNC);
    double MAF = p->params->getDoubleFlag(ARG_MAF);
    int numThreads = p->params->getIntFlag(ARG_THREAD);
    int MAX_EXTEND = ( p->params->getIntFlag(ARG_MAX_EXTEND_NSL) <= 0 ) ? nloci : p->params->getIntFlag(ARG_MAX_EXTEND_NSL);

    bool unphased = p->params->getBoolFlag(ARG_UNPHASED);

    double (*calc)(map<string, int> &, int, bool) = p->calc;

    int step = (nloci / numThreads) / (pbar->totalTicks);
    if (step == 0) step = 1;

    bool isDerived;
    string hapStr;

    for (int locus = id; locus < nloci; locus += numThreads)
    {
        if (locus % step == 0) advanceBar(*pbar, double(step));

        ihs[locus] = 0;
        //freq[locus] = MISSING;
        ihh1[locus] = MISSING;
        ihh2[locus] = MISSING;
        bool skipLocus = 0;
        //If the focal snp has MAF < MAF, then skip this locus
        if (freq[locus] < MAF || freq[locus] > 1 - MAF)
        {
            pthread_mutex_lock(&mutex_log);
            (*flog) << "WARNING: Locus " << locusName[locus]
                    << " has MAF < " << MAF << ". Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus] << "\n";
            pthread_mutex_unlock(&mutex_log);
            ihs[locus] = MISSING;
            skipLocus = 0;
            continue;
        }

        //EHH to the left of the core snp
        double current_derived_ehh = 1;
        double current_ancestral_ehh = 1;
        double previous_derived_ehh = 1;
        double previous_ancestral_ehh = 1;
        int currentLocus = locus;
        int nextLocus = locus - 1;
        double derived_ihh = 0;
        double ancestral_ihh = 0;
        
        //used for unphased analyses
        double current_notDerived_ehh = 1;
        double current_notAncestral_ehh = 1;
        double previous_notDerived_ehh = 1;
        double previous_notAncestral_ehh = 1;
        double notDerived_ihh = 0;
        double notAncestral_ihh = 0;

        //A list of all the haplotypes
        //Starts with just the focal snp and grows outward
        string *haplotypeList = new string[nhaps];
        for (int hap = 0; hap < nhaps; hap++)
        {
            //derivedCount += ( data[hap][locus] == '1' ) ? 1 : 0;
            haplotypeList[hap] = data[hap][locus];
        }

        while (current_derived_ehh > EHH_CUTOFF || current_ancestral_ehh > EHH_CUTOFF)
        {
            if (nextLocus < 0)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << EHH_CUTOFF
                        << ". ";
                if (!TRUNC){
                    skipLocus = 1;
                    (*flog) << "Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus];
                }
                (*flog) << "\n";
                pthread_mutex_unlock(&mutex_log);
                break;
            }
            else if (physicalPos[currentLocus] - physicalPos[nextLocus] > MAX_GAP)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached a gap of " << physicalPos[currentLocus] - physicalPos[nextLocus]
                        << "bp > " << MAX_GAP << "bp. Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus] << "\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }

            //Check to see if the gap between the markers is huge, if so, scale it in an ad hoc way as in
            //Voight et al. (2006)
            double scale = double(SCALE_PARAMETER) / double(physicalPos[currentLocus] - physicalPos[nextLocus]);
            if (scale > 1) scale = 1;

            currentLocus = nextLocus;
            nextLocus--;

            int numDerived = 0;
            int numAncestral = 0;

            map<string, int> ancestralHapCount;
            map<string, int> derivedHapCount;

            //used for unphased analyses
            int numHet = 0;
            map<string, int> notAncestralHapCount;
            map<string, int> notDerivedHapCount;

            for (int hap = 0; hap < nhaps; hap++)
            {
                if(unphased){
                    haplotypeList[hap] += data[hap][currentLocus];
                    hapStr = haplotypeList[hap];
                    if (data[hap][locus] == '0'){
                        //count ancestral haplotype
                        if (ancestralHapCount.count(hapStr) == 0) ancestralHapCount[hapStr] = 1;
                        else ancestralHapCount[hapStr]++;
                        numAncestral++;
                        //count non-derived hapoltype
                        if (notDerivedHapCount.count(hapStr) == 0) notDerivedHapCount[hapStr] = 1;
                        else notDerivedHapCount[hapStr]++;
                    }
                    else if (data[hap][locus] == '1'){
                        //count non-derived hapoltype
                        if (notDerivedHapCount.count(hapStr) == 0) notDerivedHapCount[hapStr] = 1;
                        else notDerivedHapCount[hapStr]++;
                        //count non-ancestral haplotype
                        if (notAncestralHapCount.count(hapStr) == 0) notAncestralHapCount[hapStr] = 1;
                        else notAncestralHapCount[hapStr]++;
                        numHet++;
                    }
                    else{//data[hap][locus] == '2'
                        //count derived hapoltype
                        if (derivedHapCount.count(hapStr) == 0) derivedHapCount[hapStr] = 1;
                        else derivedHapCount[hapStr]++;
                        numDerived++;
                        //count non-ancestral haplotype
                        if (notAncestralHapCount.count(hapStr) == 0) notAncestralHapCount[hapStr] = 1;
                        else notAncestralHapCount[hapStr]++;
                    }
                }
                else{
                    isDerived = ( data[hap][locus] == '1' ) ? 1 : 0;
                    haplotypeList[hap] += data[hap][currentLocus];
                    hapStr = haplotypeList[hap];

                    if (isDerived)
                    {
                        //count derived hapoltype
                        if (derivedHapCount.count(hapStr) == 0) derivedHapCount[hapStr] = 1;
                        else derivedHapCount[hapStr]++;
                        numDerived++;
                    }
                    else
                    {
                        //count ancestral haplotype
                        if (ancestralHapCount.count(hapStr) == 0) ancestralHapCount[hapStr] = 1;
                        else ancestralHapCount[hapStr]++;
                        numAncestral++;
                    }
                }
            }

            //We've now counted all of the unique haplotypes extending out of the core SNP
            //If locus is monomorphic, shoot a warning and skip locus
            //This probably isn't necessary any more
            if ( !unphased && (numDerived == 0 || numAncestral == 0) ) 
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: locus " << locusName[locus]
                        << " (number " << locus + 1 << ") is monomorphic. Skipping calculation at this locus.\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }

            double freqHetGT = double(numHet) / double(numDerived + numAncestral + numHet);
            //double freqAncestralGT = double(numAncestral) / double(numDerived + numAncestral + numHet);

            if ( unphased && freqHetGT > 1-MAF ) 
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: locus " << locusName[locus]
                        << " (number " << locus + 1 << ") has too many hets. Skipping calculation at this locus. "
                        << "het: " << numHet << " hom0: " << numAncestral << " hom1: " << numDerived << ".\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }

            if (current_derived_ehh > EHH_CUTOFF)
            {
                current_derived_ehh = (*calc)(derivedHapCount, numDerived, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                derived_ihh += 0.5 * scale * (current_derived_ehh + previous_derived_ehh);
                previous_derived_ehh = current_derived_ehh;
            }

            if (current_ancestral_ehh > EHH_CUTOFF)
            {
                current_ancestral_ehh = (*calc)(ancestralHapCount, numAncestral, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                ancestral_ihh += 0.5 * scale * (current_ancestral_ehh + previous_ancestral_ehh);
                previous_ancestral_ehh = current_ancestral_ehh;
            }

            if(unphased){
                current_notDerived_ehh = (*calc)(notDerivedHapCount, numAncestral + numHet, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                notDerived_ihh += 0.5 * scale * (current_notDerived_ehh + previous_notDerived_ehh);
                previous_notDerived_ehh = current_notDerived_ehh;

                current_notAncestral_ehh = (*calc)(notAncestralHapCount, numDerived + numHet, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                notAncestral_ihh += 0.5 * scale * (current_notAncestral_ehh + previous_notAncestral_ehh);
                previous_notAncestral_ehh = current_notAncestral_ehh;
            }

            //check if currentLocus is beyond MAX_EXTEND
            if (locus - currentLocus >= MAX_EXTEND) break;
        }

        delete [] haplotypeList;

        if (skipLocus == 1)
        {
            ihs[locus] = MISSING;
            skipLocus = 0;
            continue;
        }

        //calculate EHH to the right
        current_derived_ehh = 1;
        current_ancestral_ehh = 1;
        previous_derived_ehh = 1;
        previous_ancestral_ehh = 1;
        currentLocus = locus;
        nextLocus = locus + 1;
        skipLocus = 0;

        //used for unphased analyses
        current_notDerived_ehh = 1;
        current_notAncestral_ehh = 1;
        previous_notDerived_ehh = 1;
        previous_notAncestral_ehh = 1;
        notDerived_ihh = 0;
        notAncestral_ihh = 0;

        //A list of all the haplotypes
        //Starts with just the focal snp and grows outward
        haplotypeList = new string[nhaps];
        for (int hap = 0; hap < nhaps; hap++)
        {
            //char digit[2];
            //sprintf(digit, "%d", data[hap][locus]);
            haplotypeList[hap] = data[hap][locus];
        }

        while (current_ancestral_ehh > EHH_CUTOFF || current_derived_ehh > EHH_CUTOFF)
        {
            if (nextLocus > nloci - 1)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << EHH_CUTOFF
                        << ". ";
                if (!TRUNC){
                    skipLocus = 1;
                    (*flog) << "Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus];
                }
                (*flog) << "\n";
                pthread_mutex_unlock(&mutex_log);
                break;
            }
            else if (physicalPos[nextLocus] - physicalPos[currentLocus] > MAX_GAP)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached a gap of " << physicalPos[nextLocus] - physicalPos[currentLocus]
                        << "bp > " << MAX_GAP << "bp. Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus] << "\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }

            double scale = double(SCALE_PARAMETER) / double(physicalPos[nextLocus] - physicalPos[currentLocus]);
            if (scale > 1) scale = 1;

            currentLocus = nextLocus;
            nextLocus++;

            int numDerived = 0;
            int numAncestral = 0;

            map<string, int> ancestralHapCount;
            map<string, int> derivedHapCount;

            //used for unphased analyses
            int numHet = 0;
            map<string, int> notAncestralHapCount;
            map<string, int> notDerivedHapCount;

            for (int hap = 0; hap < nhaps; hap++)
            {
                if(unphased){
                    haplotypeList[hap] += data[hap][currentLocus];
                    hapStr = haplotypeList[hap];
                    if (data[hap][locus] == '0'){
                        //count ancestral haplotype
                        if (ancestralHapCount.count(hapStr) == 0) ancestralHapCount[hapStr] = 1;
                        else ancestralHapCount[hapStr]++;
                        numAncestral++;
                        //count non-derived hapoltype
                        if (notDerivedHapCount.count(hapStr) == 0) notDerivedHapCount[hapStr] = 1;
                        else notDerivedHapCount[hapStr]++;
                    }
                    else if (data[hap][locus] == '1'){
                        //count non-derived hapoltype
                        if (notDerivedHapCount.count(hapStr) == 0) notDerivedHapCount[hapStr] = 1;
                        else notDerivedHapCount[hapStr]++;
                        //count non-ancestral haplotype
                        if (notAncestralHapCount.count(hapStr) == 0) notAncestralHapCount[hapStr] = 1;
                        else notAncestralHapCount[hapStr]++;
                        numHet++;
                    }
                    else{//data[hap][locus] == '2'
                        //count derived hapoltype
                        if (derivedHapCount.count(hapStr) == 0) derivedHapCount[hapStr] = 1;
                        else derivedHapCount[hapStr]++;
                        numDerived++;
                        //count non-ancestral haplotype
                        if (notAncestralHapCount.count(hapStr) == 0) notAncestralHapCount[hapStr] = 1;
                        else notAncestralHapCount[hapStr]++;
                    }
                }
                else{
                    isDerived = ( data[hap][locus] == '1') ? 1 : 0;
                    haplotypeList[hap] += data[hap][currentLocus];
                    hapStr = haplotypeList[hap];

                    if (isDerived)
                    {
                        //count hapoltypes
                        if (derivedHapCount.count(hapStr) == 0) derivedHapCount[hapStr] = 1;
                        else derivedHapCount[hapStr]++;
                        numDerived++;
                    }
                    else //ancestral
                    {
                        if (ancestralHapCount.count(hapStr) == 0) ancestralHapCount[hapStr] = 1;
                        else ancestralHapCount[hapStr]++;
                        numAncestral++;
                    }
                }
            }

            //We've now counted all of the unique haplotypes extending out of the core SNP
            //If there are no derived alleles at a locus, shoot a warning and skip locus
            if (numDerived == 0 || numAncestral == 0)
            {
                //(*flog) << "WARNING: locus " << locusName[locus]
                //   << " (number " << locus+1 << ") is monomorphic.  Skipping calculation at this locus.\n";
                skipLocus = 1;
                break;
            }

            if (current_derived_ehh > EHH_CUTOFF)
            {
                current_derived_ehh = (*calc)(derivedHapCount, numDerived, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                derived_ihh += 0.5 * scale * (current_derived_ehh + previous_derived_ehh);
                previous_derived_ehh = current_derived_ehh;
            }

            if (current_ancestral_ehh > EHH_CUTOFF)
            {
                current_ancestral_ehh = (*calc)(ancestralHapCount, numAncestral, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                ancestral_ihh += 0.5 * scale * (current_ancestral_ehh + previous_ancestral_ehh);
                previous_ancestral_ehh = current_ancestral_ehh;
            }

            if(unphased){
                current_notDerived_ehh = (*calc)(notDerivedHapCount, numAncestral + numHet, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                notDerived_ihh += 0.5 * scale * (current_notDerived_ehh + previous_notDerived_ehh);
                previous_notDerived_ehh = current_notDerived_ehh;

                current_notAncestral_ehh = (*calc)(notAncestralHapCount, numDerived + numHet, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                notAncestral_ihh += 0.5 * scale * (current_notAncestral_ehh + previous_notAncestral_ehh);
                previous_notAncestral_ehh = current_notAncestral_ehh;
            }

            //check if currentLocus is beyond MAX_EXTEND
            if (currentLocus - locus >= MAX_EXTEND) break;

        }

        delete [] haplotypeList;

        if (skipLocus == 1)
        {
            ihs[locus] = MISSING;
            skipLocus = 0;
            continue;
        }

        if (ihs[locus] != MISSING)
        {
            if(unphased){
                ihh1[locus] = log10(derived_ihh / notDerived_ihh);
                ihh2[locus] = log10(ancestral_ihh / notAncestral_ihh);
                ihs[locus] = (ihh1[locus] > ihh2[locus]) ? ihh1[locus] : 0-ihh2[locus];
            }
            else{
                ihh1[locus] = derived_ihh;
                ihh2[locus] = ancestral_ihh;
                ihs[locus] = log10(derived_ihh / ancestral_ihh);
                //freq[locus] = double(derivedCount) / double(nhaps);
            }
        }
    }
}


void calc_soft_ihs(void *order)
{
    work_order_t *p = (work_order_t *)order;
    char **data = p->hapData->data;
    int nloci = p->hapData->nloci;
    int nhaps = p->hapData->nhaps;
    int *physicalPos = p->mapData->physicalPos;
    double *geneticPos = p->mapData->geneticPos;
    string *locusName = p->mapData->locusName;
    //int start = p->first_index;
    //int stop = p->last_index;
    int id = p->id;
    double *h1 = p->ihh1;
    double *h2dh1 = p->ihh2;
    double *h12 = p->ihs;
    //double *freq = p->freq;
    ofstream *flog = p->flog;
    Bar *pbar = p->bar;

    int SCALE_PARAMETER = p->params->getIntFlag(ARG_GAP_SCALE);
    int MAX_GAP = p->params->getIntFlag(ARG_MAX_GAP);
    double EHH_CUTOFF = p->params->getDoubleFlag(ARG_CUTOFF);
    //bool ALT = p->params->getBoolFlag(ARG_ALT);
    //double MAF = p->params->getDoubleFlag(ARG_MAF);
    //bool WAGH = p->params->getBoolFlag(ARG_WAGH);
    int numThreads = p->params->getIntFlag(ARG_THREAD);
    bool TRUNC = p->params->getBoolFlag(ARG_TRUNC);
    int MAX_EXTEND = ( p->params->getIntFlag(ARG_MAX_EXTEND) <= 0 ) ? physicalPos[nloci - 1] - physicalPos[0] : p->params->getIntFlag(ARG_MAX_EXTEND);

    int step = (nloci / numThreads) / (pbar->totalTicks);
    if (step == 0) step = 1;

    for (int locus = id; locus < nloci; locus += numThreads)
    {
        if (locus % step == 0) advanceBar(*pbar, double(step));

        //freq[locus] = MISSING;
        h1[locus] = -1;
        h2dh1[locus] = -1;
        h12[locus] = -1;

        //EHH to the left of the core snp
        double current_ehh1 = 1;
        double current_ehh2d1 = 1;
        double current_ehh12 = 1;

        double previous_ehh1 = 1;
        double previous_ehh2d1 = 1;
        double previous_ehh12 = 1;

        int currentLocus = locus;
        int nextLocus = locus - 1;
        bool skipLocus = 0;
        double ihh1 = 0;
        double ihh2d1 = 0;
        double ihh12 = 0;
        //double derivedCount = 0;
        //A list of all the haplotypes
        //Starts with just the focal snp and grows outward
        map<string, int> tempHapCount;
        string *haplotypeList = new string[nhaps];
        for (int hap = 0; hap < nhaps; hap++)
        {
            //derivedCount += data[hap][locus];
            haplotypeList[hap] = data[hap][locus];
            string hapStr = haplotypeList[hap];
            //count hapoltype freqs
            if (tempHapCount.count(hapStr) == 0) tempHapCount[hapStr] = 1;
            else tempHapCount[hapStr]++;
        }

        triplet_t res = calculateSoft(tempHapCount, nhaps);
        current_ehh1 = res.h1;
        current_ehh2d1 = res.h2dh1;
        current_ehh12 = res.h12;
        previous_ehh1 = res.h1;
        previous_ehh2d1 = res.h2dh1;
        previous_ehh12 = res.h12;

        while (current_ehh1 > EHH_CUTOFF)
        {
            if (nextLocus < 0)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << EHH_CUTOFF << ". ";
                if (!TRUNC){
                    skipLocus = 1;
                    (*flog) << "Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus];
                }
                (*flog) << "\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }
            else if (physicalPos[currentLocus] - physicalPos[nextLocus] > MAX_GAP)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached a gap of " << physicalPos[currentLocus] - physicalPos[nextLocus] << "bp > " << MAX_GAP 
                << "bp. Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus] << "\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }

            //Check to see if the gap between the markers is huge, if so, scale it in an ad hoc way as in
            //Voight et al. (2006)
            double scale = double(SCALE_PARAMETER) / double(physicalPos[currentLocus] - physicalPos[nextLocus]);
            if (scale > 1) scale = 1;

            currentLocus = nextLocus;
            nextLocus--;

            map<string, int> hapCount;

            for (int hap = 0; hap < nhaps; hap++)
            {
                //build haplotype string
                haplotypeList[hap] += data[hap][currentLocus];
                string hapStr = haplotypeList[hap];

                //count hapoltype freqs
                if (hapCount.count(hapStr) == 0) hapCount[hapStr] = 1;
                else hapCount[hapStr]++;
            }

            //We've now counted all of the unique haplotypes extending out of the core SNP
            res = calculateSoft(hapCount, nhaps);
            current_ehh1 = res.h1;
            current_ehh2d1 = res.h2dh1;
            current_ehh12 = res.h12;

            //directly calculate integral, iteratively
            //Trapezoid rule
            ihh1 += 0.5 * scale * (geneticPos[currentLocus + 1] - geneticPos[currentLocus]) * (current_ehh1 + previous_ehh1);
            previous_ehh1 = current_ehh1;
            ihh2d1 += 0.5 * scale * (geneticPos[currentLocus + 1] - geneticPos[currentLocus]) * (current_ehh2d1 + previous_ehh2d1);
            previous_ehh2d1 = current_ehh2d1;
            ihh12 += 0.5 * scale * (geneticPos[currentLocus + 1] - geneticPos[currentLocus]) * (current_ehh12 + previous_ehh12);
            previous_ehh12 = current_ehh12;

            //check if currentLocus is beyond 1Mb
            if (physicalPos[locus] - physicalPos[currentLocus] >= MAX_EXTEND) break;
        }

        delete [] haplotypeList;

        if (skipLocus == 1)
        {
            h1[locus] = MISSING;
            h2dh1[locus] = MISSING;
            h12[locus] = MISSING;
            skipLocus = 0;
            continue;
        }

        //calculate EHH to the right
        current_ehh1 = 1;
        current_ehh2d1 = 1;
        current_ehh12 = 1;

        previous_ehh1 = 1;
        previous_ehh2d1 = 1;
        previous_ehh12 = 1;

        ihh1 = 0;
        ihh2d1 = 0;
        ihh12 = 0;

        currentLocus = locus;
        nextLocus = locus + 1;
        skipLocus = 0;
        //A list of all the haplotypes
        //Starts with just the focal snp and grows outward
        tempHapCount.clear();
        haplotypeList = new string[nhaps];
        for (int hap = 0; hap < nhaps; hap++)
        {
            haplotypeList[hap] = data[hap][locus];
            string hapStr = haplotypeList[hap];
            //count hapoltype freqs
            if (tempHapCount.count(hapStr) == 0) tempHapCount[hapStr] = 1;
            else tempHapCount[hapStr]++;
        }

        res = calculateSoft(tempHapCount, nhaps);
        current_ehh1 = res.h1;
        current_ehh2d1 = res.h2dh1;
        current_ehh12 = res.h12;
        previous_ehh1 = res.h1;
        previous_ehh2d1 = res.h2dh1;
        previous_ehh12 = res.h12;

        while (current_ehh1 > EHH_CUTOFF)
        {
            if (nextLocus > nloci - 1)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << EHH_CUTOFF << ". ";
                if (!TRUNC){
                    skipLocus = 1;
                    (*flog) << "Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus];
                }
                (*flog) << "\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }
            else if (physicalPos[nextLocus] - physicalPos[currentLocus] > MAX_GAP)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached a gap of " << physicalPos[nextLocus] - physicalPos[currentLocus] << "bp > " << MAX_GAP << "bp. Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus] << "\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }

            double scale = double(SCALE_PARAMETER) / double(physicalPos[nextLocus] - physicalPos[currentLocus]);
            if (scale > 1) scale = 1;

            currentLocus = nextLocus;
            nextLocus++;

            map<string, int> hapCount;

            for (int hap = 0; hap < nhaps; hap++)
            {
                //build haplotype string
                haplotypeList[hap] += data[hap][currentLocus];
                string hapStr = haplotypeList[hap];

                //count hapoltype freqs
                if (hapCount.count(hapStr) == 0) hapCount[hapStr] = 1;
                else hapCount[hapStr]++;
            }

            //We've now counted all of the unique haplotypes extending out of the core SNP
            res = calculateSoft(hapCount, nhaps);
            current_ehh1 = res.h1;
            current_ehh2d1 = res.h2dh1;
            current_ehh12 = res.h12;

            //directly calculate integral, iteratively
            //Trapezoid rule
            ihh1 += 0.5 * scale * (geneticPos[currentLocus] - geneticPos[currentLocus - 1]) * (current_ehh1 + previous_ehh1);
            previous_ehh1 = current_ehh1;
            ihh2d1 += 0.5 * scale * (geneticPos[currentLocus] - geneticPos[currentLocus - 1]) * (current_ehh2d1 + previous_ehh2d1);
            previous_ehh2d1 = current_ehh2d1;
            ihh12 += 0.5 * scale * (geneticPos[currentLocus] - geneticPos[currentLocus - 1]) * (current_ehh12 + previous_ehh12);
            previous_ehh12 = current_ehh12;

            //check if currentLocus is beyond 1Mb
            if (physicalPos[currentLocus] - physicalPos[locus] >= MAX_EXTEND) break;
        }

        delete [] haplotypeList;

        if (skipLocus == 1)
        {
            h1[locus] = MISSING;
            h2dh1[locus] = MISSING;
            h12[locus] = MISSING;
            skipLocus = 0;
            continue;
        }

        if (h12[locus] != MISSING)
        {
            h1[locus] = ihh1;
            h2dh1[locus] = ihh2d1;
            h12[locus] = ihh12;
            //freq[locus] = double(derivedCount) / double(nhaps);
        }
    }
}


void calc_xpihh(void *order)
{
    work_order_t *p = (work_order_t *)order;

    char **data1 = p->hapData1->data;
    int nhaps1 = p->hapData1->nhaps;
    double *ihh1 = p->ihh1;
    //double *freq1 = p->freq1;

    char **data2 = p->hapData2->data;
    int nhaps2 = p->hapData2->nhaps;
    double *ihh2 = p->ihh2;
    //double *freq2 = p->freq2;

    int nloci = p->mapData->nloci;
    int *physicalPos = p->mapData->physicalPos;
    double *geneticPos = p->mapData->geneticPos;
    string *locusName = p->mapData->locusName;

    int id = p->id;

    ofstream *flog = p->flog;
    Bar *pbar = p->bar;

    bool TRUNC = p->params->getBoolFlag(ARG_TRUNC);
    int SCALE_PARAMETER = p->params->getIntFlag(ARG_GAP_SCALE);
    int MAX_GAP = p->params->getIntFlag(ARG_MAX_GAP);
    double EHH_CUTOFF = p->params->getDoubleFlag(ARG_CUTOFF);
    bool ALT = p->params->getBoolFlag(ARG_ALT);
    bool WAGH = p->params->getBoolFlag(ARG_WAGH);
    int numThreads = p->params->getIntFlag(ARG_THREAD);
    bool CALC_XPNSL = p->params->getBoolFlag(ARG_XPNSL);

    int MAX_EXTEND;
    if (!CALC_XPNSL){
        MAX_EXTEND = ( p->params->getIntFlag(ARG_MAX_EXTEND) <= 0 ) ? physicalPos[nloci - 1] - physicalPos[0] : p->params->getIntFlag(ARG_MAX_EXTEND);
    }
    else{
        MAX_EXTEND = ( p->params->getIntFlag(ARG_MAX_EXTEND_NSL) <= 0 ) ? geneticPos[nloci - 1] - geneticPos[0] : p->params->getIntFlag(ARG_MAX_EXTEND_NSL);
    }
    int step = (nloci / numThreads) / (pbar->totalTicks);
    if (step == 0) step = 1;

    string hapStr;

    for (int locus = id; locus < nloci; locus += numThreads)
    {
        if (locus % step == 0) advanceBar(*pbar, double(step));

        ihh1[locus] = 0;
        ihh2[locus] = 0;

        //freq1[locus] = MISSING;
        //freq2[locus] = MISSING;

        //EHH to the left of the core snp
        double current_pooled_ehh = 1;
        double previous_pooled_ehh = 1;
        double derivedCountPooled = 0;

        double current_pop1_ehh = 1;
        double previous_pop1_ehh = 1;
        double ihhPop1 = 0;
        double derivedCount1 = 0;

        double current_pop2_ehh = 1;
        double previous_pop2_ehh = 1;
        double ihhPop2 = 0;
        double derivedCount2 = 0;

        int currentLocus = locus;
        int nextLocus = locus - 1;
        bool skipLocus = 0;

        //A list of all the haplotypes
        //Starts with just the focal snp and grows outward
        string *haplotypeList1, *haplotypeList2, *haplotypeListPooled;
        haplotypeList1 = new string[nhaps1];
        haplotypeList2 = new string[nhaps2];
        haplotypeListPooled = new string[nhaps1 + nhaps2];
        for (int hap = 0; hap < nhaps1 + nhaps2; hap++)
        {
            //char digit[2];
            //Pop1
            if (hap < nhaps1)
            {
                //sprintf(digit, "%d", data1[hap][locus]);
                haplotypeList1[hap] = data1[hap][locus];
                haplotypeListPooled[hap] = data1[hap][locus];
                derivedCount1 += ( data1[hap][locus] == '1') ? 1 : 0;
            }
            //Pop2
            else
            {
                //sprintf(digit, "%d", data2[hap - nhaps1][locus]);
                haplotypeList2[hap - nhaps1] = data2[hap - nhaps1][locus];
                haplotypeListPooled[hap] = data2[hap - nhaps1][locus];
                derivedCount2 += ( data2[hap - nhaps1][locus] == '1' ) ? 1 : 0;
            }
        }

        derivedCountPooled = derivedCount1 + derivedCount2;

        //when calculating xp-ehh, ehh does not necessarily start at 1
        if (ALT)
        {
            double f = double(derivedCount1) / double(nhaps1);
            current_pop1_ehh = f * f + (1 - f) * (1 - f);
            previous_pop1_ehh = current_pop1_ehh;

            f = double(derivedCount2) / double(nhaps2);
            current_pop2_ehh = f * f + (1 - f) * (1 - f);
            previous_pop2_ehh = current_pop2_ehh;

            f = double(derivedCountPooled) / double(nhaps1 + nhaps2);
            current_pooled_ehh = f * f + (1 - f) * (1 - f);
            previous_pooled_ehh = current_pooled_ehh;
        }
        else
        {
            if (WAGH)
            {

                current_pop1_ehh = (derivedCount1 > 1) ? nCk(derivedCount1,2) / (nCk(derivedCount1,2)+nCk(nhaps1-derivedCount1,2)) : 0;
                current_pop1_ehh += (nhaps1 - derivedCount1 > 1) ? nCk(nhaps1-derivedCount1,2) / (nCk(derivedCount1,2)+nCk(nhaps1-derivedCount1,2)) : 0;
                previous_pop1_ehh = current_pop1_ehh;

                current_pop2_ehh = (derivedCount2 > 1) ? nCk(derivedCount2, 2) / (nCk(derivedCount2,2)+nCk(nhaps2-derivedCount2,2)) : 0;
                current_pop2_ehh += (nhaps2 - derivedCount2 > 1) ? nCk(nhaps2 - derivedCount2, 2) / (nCk(derivedCount2,2)+nCk(nhaps2-derivedCount2,2)) : 0;
                previous_pop2_ehh = current_pop2_ehh;

            }
            else
            {
                current_pop1_ehh = (derivedCount1 > 1) ? nCk(derivedCount1, 2) / nCk(nhaps1, 2) : 0;
                current_pop1_ehh += (nhaps1 - derivedCount1 > 1) ? nCk(nhaps1 - derivedCount1, 2) / nCk(nhaps1, 2) : 0;
                previous_pop1_ehh = current_pop1_ehh;

                current_pop2_ehh = (derivedCount2 > 1) ? nCk(derivedCount2, 2) / nCk(nhaps2, 2) : 0;
                current_pop2_ehh += (nhaps2 - derivedCount2 > 1) ? nCk(nhaps2 - derivedCount2, 2) / nCk(nhaps2, 2) : 0;
                previous_pop2_ehh = current_pop2_ehh;        

            }

            current_pooled_ehh = (derivedCountPooled > 1) ? nCk(derivedCountPooled, 2) / nCk(nhaps1 + nhaps2, 2) : 0;
            current_pooled_ehh += (nhaps1 + nhaps2 - derivedCountPooled > 1) ? nCk(nhaps1 + nhaps2 - derivedCountPooled, 2) / nCk(nhaps1 + nhaps2, 2) : 0;
            previous_pooled_ehh = current_pooled_ehh;    
        }

        while (current_pooled_ehh > EHH_CUTOFF)
        {
            if (nextLocus < 0)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << EHH_CUTOFF
                        << ". ";
                if (!TRUNC){
                    skipLocus = 1;
                    (*flog) << "Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus];
                }
                (*flog) << "\n";
                pthread_mutex_unlock(&mutex_log);
                break;
            }
            else if (physicalPos[currentLocus] - physicalPos[nextLocus] > MAX_GAP)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached a gap of " << physicalPos[currentLocus] - physicalPos[nextLocus] << "bp > " << MAX_GAP 
                << "bp. Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus] << "\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }

            //Check to see if the gap between the markers is huge, if so, scale it in an ad hoc way as in
            //Voight, et al. paper
            double scale = double(SCALE_PARAMETER) / double(physicalPos[currentLocus] - physicalPos[nextLocus]);
            if (scale > 1) scale = 1;

            currentLocus = nextLocus;
            nextLocus--;

            map<string, int> hapCount1;
            map<string, int> hapCount2;
            map<string, int> hapCountPooled;

            //build haplotype strings
            for (int hap = 0; hap < nhaps1 + nhaps2; hap++)
            {
                //char digit[2];
                if (hap < nhaps1) //Pop1
                {
                    //sprintf(digit, "%d", data1[hap][currentLocus]);
                    haplotypeList1[hap] += data1[hap][currentLocus];
                    hapStr = haplotypeList1[hap];
                    if (hapCount1.count(hapStr) == 0) hapCount1[hapStr] = 1;
                    else hapCount1[hapStr]++;

                    //Pooled
                    haplotypeListPooled[hap] += data1[hap][currentLocus];
                    hapStr = haplotypeListPooled[hap];
                    if (hapCountPooled.count(hapStr) == 0) hapCountPooled[hapStr] = 1;
                    else hapCountPooled[hapStr]++;
                }
                else //Pop2
                {
                    //sprintf(digit, "%d", data2[hap - nhaps1][currentLocus]);
                    haplotypeList2[hap - nhaps1] += data2[hap - nhaps1][currentLocus];
                    hapStr = haplotypeList2[hap - nhaps1];
                    if (hapCount2.count(hapStr) == 0) hapCount2[hapStr] = 1;
                    else hapCount2[hapStr]++;

                    //Pooled
                    haplotypeListPooled[hap] += data2[hap - nhaps1][currentLocus];
                    hapStr = haplotypeListPooled[hap];
                    if (hapCountPooled.count(hapStr) == 0) hapCountPooled[hapStr] = 1;
                    else hapCountPooled[hapStr]++;
                }
            }

            if (WAGH)
            {
                current_pop1_ehh = calculateHomozygosity_Wagh(hapCount1,nhaps1,derivedCount1);
                current_pop2_ehh = calculateHomozygosity_Wagh(hapCount2,nhaps2,derivedCount2);

            }
            else
            {
                current_pop1_ehh = calculateHomozygosity(hapCount1, nhaps1, ALT);
                current_pop2_ehh = calculateHomozygosity(hapCount2, nhaps2, ALT);

            }

                       
            current_pooled_ehh = calculateHomozygosity(hapCountPooled, nhaps1 + nhaps2, ALT);

            //directly calculate ihh, iteratively
            //Trapezoid rule
            ihhPop1 += 0.5 * scale * (geneticPos[currentLocus + 1] - geneticPos[currentLocus]) * (current_pop1_ehh + previous_pop1_ehh);
            previous_pop1_ehh = current_pop1_ehh;

            ihhPop2 += 0.5 * scale * (geneticPos[currentLocus + 1] - geneticPos[currentLocus]) * (current_pop2_ehh + previous_pop2_ehh);
            previous_pop2_ehh = current_pop2_ehh;

            previous_pooled_ehh = current_pooled_ehh;

            //check if currentLocus is beyond MAX_EXTEND
            if (!CALC_XPNSL && physicalPos[locus] - physicalPos[currentLocus] >= MAX_EXTEND) break;
            if (CALC_XPNSL && geneticPos[locus] - geneticPos[currentLocus] >= MAX_EXTEND) break;
        }

        delete [] haplotypeList1;
        delete [] haplotypeList2;
        delete [] haplotypeListPooled;

        if (skipLocus == 1)
        {
            ihh1[locus] = MISSING;
            ihh2[locus] = MISSING;
            skipLocus = 0;
            continue;
        }

        //calculate EHH to the right

        current_pooled_ehh = 1;
        previous_pooled_ehh = 1;
        derivedCountPooled = 0;

        current_pop1_ehh = 1;
        previous_pop1_ehh = 1;
        derivedCount1 = 0;

        current_pop2_ehh = 1;
        previous_pop2_ehh = 1;
        derivedCount2 = 0;

        currentLocus = locus;
        nextLocus = locus + 1;
        skipLocus = 0;

        //A list of all the haplotypes
        //Starts with just the focal snp and grows outward
        haplotypeList1 = new string[nhaps1];
        haplotypeList2 = new string[nhaps2];
        haplotypeListPooled = new string[nhaps1 + nhaps2];
        for (int hap = 0; hap < nhaps1 + nhaps2; hap++)
        {
            //Pop1
            if (hap < nhaps1)
            {
                //sprintf(digit, "%d", data1[hap][locus]);
                haplotypeList1[hap] = data1[hap][locus];
                derivedCount1 += ( data1[hap][locus] == '1' ) ? 1 : 0;

                //Pooled
                haplotypeListPooled[hap] = data1[hap][locus];
            }
            //Pop2
            else
            {
                //sprintf(digit, "%d", data2[hap - nhaps1][locus]);
                haplotypeList2[hap - nhaps1] = data2[hap - nhaps1][locus];
                derivedCount2 += ( data2[hap - nhaps1][locus] == '1' ) ? 1 : 0;

                //Pooled
                haplotypeListPooled[hap] = data2[hap - nhaps1][locus];
            }
        }

        derivedCountPooled = derivedCount1 + derivedCount2;

        //when calculating xp-ehh, ehh does not necessarily start at 1
        if (ALT)
        {
            double f = double(derivedCount1) / double(nhaps1);
            current_pop1_ehh = f * f + (1 - f) * (1 - f);
            previous_pop1_ehh = current_pop1_ehh;

            f = double(derivedCount2) / double(nhaps2);
            current_pop2_ehh = f * f + (1 - f) * (1 - f);
            previous_pop2_ehh = current_pop2_ehh;

            f = double(derivedCountPooled) / double(nhaps1 + nhaps2);
            current_pooled_ehh = f * f + (1 - f) * (1 - f);
            previous_pooled_ehh = current_pooled_ehh;
        }
        else
        {

            if (WAGH)
            {
                current_pop1_ehh = (derivedCount1 > 1) ? nCk(derivedCount1,2) / (nCk(derivedCount1,2)+nCk(nhaps1-derivedCount1,2)) : 0;
                current_pop1_ehh += (nhaps1 - derivedCount1 > 1) ? nCk(nhaps1-derivedCount1,2) / (nCk(derivedCount1,2)+nCk(nhaps1-derivedCount1,2)) : 0;
                previous_pop1_ehh = current_pop1_ehh;

                current_pop2_ehh = (derivedCount2 > 1) ? nCk(derivedCount2, 2) / (nCk(derivedCount2,2)+nCk(nhaps2-derivedCount2,2)) : 0;
                current_pop2_ehh += (nhaps2 - derivedCount2 > 1) ? nCk(nhaps2 - derivedCount2, 2) / (nCk(derivedCount2,2)+nCk(nhaps2-derivedCount2,2)) : 0;
                previous_pop2_ehh = current_pop2_ehh;                  
            }
            else
            {
                current_pop1_ehh = (derivedCount1 > 1) ? nCk(derivedCount1, 2) / nCk(nhaps1, 2) : 0;
                current_pop1_ehh += (nhaps1 - derivedCount1 > 1) ? nCk(nhaps1 - derivedCount1, 2) / nCk(nhaps1, 2) : 0;
                previous_pop1_ehh = current_pop1_ehh;

                current_pop2_ehh = (derivedCount2 > 1) ? nCk(derivedCount2, 2) / nCk(nhaps2, 2) : 0;
                current_pop2_ehh += (nhaps2 - derivedCount2 > 1) ? nCk(nhaps2 - derivedCount2, 2) / nCk(nhaps2, 2) : 0; 
                previous_pop2_ehh = current_pop2_ehh;                
            }
            // current_pop1_ehh = (derivedCount1 > 1) ? nCk(derivedCount1, 2) / nCk(nhaps1, 2) : 0;
            // current_pop1_ehh += (nhaps1 - derivedCount1 > 1) ? nCk(nhaps1 - derivedCount1, 2) / nCk(nhaps1, 2) : 0;
            // previous_pop1_ehh = current_pop1_ehh;

            // current_pop2_ehh = (derivedCount2 > 1) ? nCk(derivedCount2, 2) / nCk(nhaps2, 2) : 0;
            // current_pop2_ehh += (nhaps2 - derivedCount2 > 1) ? nCk(nhaps2 - derivedCount2, 2) / nCk(nhaps2, 2) : 0;
            // previous_pop2_ehh = current_pop2_ehh;

            // current_pop1_ehh = (derivedCount1 > 1) ? nCk(derivedCount1, 2) / nCk(nhaps1, 2) : 0;
            // current_pop1_ehh += (nhaps1 - derivedCount1 > 1) ? nCk(nhaps1 - derivedCount1, 2) / nCk(nhaps1, 2) : 0;
            // current_pop1_ehh = (derivedCount1 > 1) ? nCk(derivedCount1,2) / (nCk(derivedCount1,2)+nCk(nhaps1-derivedCount1,2)) : 0;
            // current_pop1_ehh += (nhaps1 - derivedCount1 > 1) ? nCk(nhaps1-derivedCount1,2) / (nCk(derivedCount1,2)+nCk(nhaps1-derivedCount1,2)) : 0;
            // previous_pop1_ehh = current_pop1_ehh;

     

            // current_pop2_ehh = (derivedCount2 > 1) ? nCk(derivedCount2, 2) / nCk(nhaps2, 2) : 0;
            // current_pop2_ehh += (nhaps2 - derivedCount2 > 1) ? nCk(nhaps2 - derivedCount2, 2) / nCk(nhaps2, 2) : 0;
            // current_pop2_ehh = (derivedCount2 > 1) ? nCk(derivedCount2, 2) / (nCk(derivedCount2,2)+nCk(nhaps2-derivedCount2,2)) : 0;
            // current_pop2_ehh += (nhaps2 - derivedCount2 > 1) ? nCk(nhaps2 - derivedCount2, 2) / (nCk(derivedCount2,2)+nCk(nhaps2-derivedCount2,2)) : 0;
            // previous_pop2_ehh = current_pop2_ehh;            

            current_pooled_ehh = (derivedCountPooled > 1) ? nCk(derivedCountPooled, 2) / nCk(nhaps1 + nhaps2, 2) : 0;
            current_pooled_ehh += (nhaps1 + nhaps2 - derivedCountPooled > 1) ? nCk(nhaps1 + nhaps2 - derivedCountPooled, 2) / nCk(nhaps1 + nhaps2, 2) : 0;
            previous_pooled_ehh = current_pooled_ehh;
        }


        while (current_pooled_ehh > EHH_CUTOFF)
        {
            if (nextLocus > nloci - 1)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << EHH_CUTOFF
                        << ". ";
                if (!TRUNC){
                    skipLocus = 1;
                    (*flog) << "Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus];
                }
                (*flog) << "\n";
                pthread_mutex_unlock(&mutex_log);
                break;
            }
            else if (physicalPos[nextLocus] - physicalPos[currentLocus] > MAX_GAP)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached a gap of " << physicalPos[nextLocus] - physicalPos[currentLocus] << "bp > " << MAX_GAP << "bp. Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus] << "\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }

            double scale = double(SCALE_PARAMETER) / double(physicalPos[nextLocus] - physicalPos[currentLocus]);
            if (scale > 1) scale = 1;

            currentLocus = nextLocus;
            nextLocus++;

            map<string, int> hapCount1;
            map<string, int> hapCount2;
            map<string, int> hapCountPooled;

            //build haplotype strings
            for (int hap = 0; hap < nhaps1 + nhaps2; hap++)
            {
                //char digit[2];
                //string hapStr;

                //pop1
                if (hap < nhaps1)
                {
                    haplotypeList1[hap] += data1[hap][currentLocus];
                    hapStr = haplotypeList1[hap];
                    if (hapCount1.count(hapStr) == 0) hapCount1[hapStr] = 1;
                    else hapCount1[hapStr]++;

                    //Pooled
                    haplotypeListPooled[hap] += data1[hap][currentLocus];
                    hapStr = haplotypeListPooled[hap];
                    if (hapCountPooled.count(hapStr) == 0) hapCountPooled[hapStr] = 1;
                    else hapCountPooled[hapStr]++;
                }
                //Pop2
                else
                {
                    haplotypeList2[hap - nhaps1] += data2[hap - nhaps1][currentLocus];
                    hapStr = haplotypeList2[hap - nhaps1];
                    if (hapCount2.count(hapStr) == 0) hapCount2[hapStr] = 1;
                    else hapCount2[hapStr]++;

                    //Pooled
                    haplotypeListPooled[hap] += data2[hap - nhaps1][currentLocus];
                    hapStr = haplotypeListPooled[hap];
                    if (hapCountPooled.count(hapStr) == 0) hapCountPooled[hapStr] = 1;
                    else hapCountPooled[hapStr]++;
                }
            }

            if (WAGH)
            {
                current_pop1_ehh = calculateHomozygosity_Wagh(hapCount1,nhaps1,derivedCount1);
                current_pop2_ehh = calculateHomozygosity_Wagh(hapCount2,nhaps2,derivedCount2);

            }
            else
            {

                current_pop1_ehh = calculateHomozygosity(hapCount1, nhaps1, ALT);
                current_pop2_ehh = calculateHomozygosity(hapCount2, nhaps2, ALT);
            }


            current_pooled_ehh = calculateHomozygosity(hapCountPooled, nhaps1 + nhaps2, ALT);

            //directly calculate ihh1, iteratively
            //Trapezoid rule
            ihhPop1 += 0.5 * scale * (geneticPos[currentLocus] - geneticPos[currentLocus - 1]) * (current_pop1_ehh + previous_pop1_ehh);
            previous_pop1_ehh = current_pop1_ehh;

            ihhPop2 += 0.5 * scale * (geneticPos[currentLocus] - geneticPos[currentLocus - 1]) * (current_pop2_ehh + previous_pop2_ehh);
            previous_pop2_ehh = current_pop2_ehh;

            previous_pooled_ehh = current_pooled_ehh;

            //check if currentLocus is beyond 1Mb
            if (!CALC_XPNSL && physicalPos[currentLocus] - physicalPos[locus] >= MAX_EXTEND) break;
            if (CALC_XPNSL && geneticPos[currentLocus] - geneticPos[locus] >= MAX_EXTEND) break;
        }

        delete [] haplotypeList1;
        delete [] haplotypeList2;
        delete [] haplotypeListPooled;

        if (skipLocus == 1)
        {
            ihh1[locus] = MISSING;
            ihh2[locus] = MISSING;
            skipLocus = 0;
            continue;
        }

        if (ihh1[locus] != MISSING)
        {
            ihh1[locus] = ihhPop1;
            //freq1[locus] = double(derivedCount1) / double(nhaps1);
        }

        if (ihh2[locus] != MISSING)
        {
            ihh2[locus] = ihhPop2;
            //freq2[locus] = double(derivedCount2) / double(nhaps2);
        }
    }
}

double calculateHomozygosity_Wagh(map<string, int> &count, int total, int derivedCount)
{
    double freq = 0;
    double homozygosity = 0;
    map<string, int>::iterator it;
    for (it = count.begin(); it != count.end(); it++)
    {
        homozygosity += (it->second > 1) ? nCk(it->second,2)/(nCk(derivedCount, 2) + nCk(total-derivedCount, 2)) : 0;
    }
        
    return homozygosity;
}


double calculateHomozygosity(map<string, int> &count, int total, bool ALT) // Called by XP-EHH
{
    double freq = 0;
    double homozygosity = 0;
    map<string, int>::iterator it;
    for (it = count.begin(); it != count.end(); it++)
    {
        if (ALT)
        {
            freq = double(it->second) / double(total);
            homozygosity += freq * freq;
        }

        else
        {
            homozygosity += (it->second > 1) ? nCk(it->second, 2) / nCk(total, 2) : 0;
        }
    }

    return homozygosity;
}

//Have to modify the function to handle the arbitrary declaration
//of EHH1K_VALUES
//If we pass, the MAX of that vector, plus the vector iteself
//We should be able to track up to the max in an array (up to what we do now for k = 2)
//return an array that is pre-allocated to the side of EKK1k_VALUES,
//and the index of that array corresponds to the indicies of EHH1K_VALUES
triplet_t calculateSoft(map<string, int> &count, int total)
{
    triplet_t res;
    res.h1 = 0;
    res.h2dh1 = 0;
    res.h12 = 0;

    double first = 0;
    double second = 0;
    //double freq = 0;
    double homozygosity = 0;
    map<string, int>::iterator it;
    for (it = count.begin(); it != count.end(); it++)
    {
        homozygosity += (it->second > 1) ? nCk(it->second, 2) / nCk(total, 2) : 0;
        //We can either do what is effectively a bubble sort here for the first K
        //sorted values
        //Or we can do a complete sort outside of the for loop
        //Will have to think about the scaling issues here
        if (it->second > first)
        {
            second = first;
            first = it->second;
        }
        else if (it->second > second)
        {
            second = it->second;
        }
    }

    //here we have a forloop over the EHH1k_VALUES elements
    double firstFreq = (first > 1) ? nCk(first, 2) / nCk(total, 2) : 0;
    double secondFreq = (second > 1) ? nCk(second, 2) / nCk(total, 2) : 0;
    double comboFreq = ((first + second) > 1) ? nCk((first + second), 2) / nCk(total, 2) : 0;

    res.h1 = homozygosity;
    res.h2dh1 = (homozygosity - firstFreq) / homozygosity;
    res.h12 = homozygosity - firstFreq - secondFreq + comboFreq;

    return res;
}

