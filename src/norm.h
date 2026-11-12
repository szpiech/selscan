#ifndef __NORM_H__
#define __NORM_H__


/* norm -- a program for downstream analysis of iHS scores calculated by selscan
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

#include "version.h"
#include <iostream>
#include <fstream>
#include <map>
#include <cmath>
#include <cstdlib>
#include <cctype>
#include <vector>
#include <cstring>
#include <sstream>
#include "param_t.h"


//
#include <iostream>
#include <vector>
#include <algorithm> // for std::shuffle
#include <numeric>   // for std::accumulate
#include <ctime>     // for seeding

#include "gene.h"
#include "selscan-numeric.h"

using namespace std;

// added afterwards

struct BinData {
    std::vector<double> bins;
    std::vector<double> means;
    std::vector<double> variances;
};


class SelscanNorm{
    public:

    GeneAnalyzer genex;
    int WINOFFSET;

    
    const string VERSION = SELSCAN_VERSION;   // was "3.0.0" here while selscan said 3.0.1

    const string PREAMBLE = " -- norm: a program for downstream analysis of selscan output\n\
    Source code and binaries can be found at\n\
    \t<https://www.github.com/szpiech/selscan>\n\
    \n\
    Citations:\n\
    \n\
    selscan: ZA Szpiech and RD Hernandez (2014) MBE 31: 2824-2827.\n\
            ZA Szpiech (2024) Bioinformatics 40: btae006.\n\
            A Rahman, TQ Smith, and ZA Szpiech (2025). MBE, 42(11), msaf275.\n\
            A Rahman, TQ Smith, and ZA Szpiech (2026). HPGG, 6(2), 0005. https://doi.org/10.47248/hpgg2606020005.\n\
    iHH12: R Torres et al. (2018) PLoS Genetics 15: e1007898.\n\
        N Garud et al. (2015) PLoS Genetics 11: 1–32.\n\
    nSL: A Ferrer-Admetlla et al. (2014) MBE 31: 1275-1291.\n\
    XP-nSL: Szpiech et al. (2021) Evol Lett 5: 408-421.\n\
    XP-EHH: PC Sabeti et al. (2007) Nature 449: 913–918.\n\
            K Wagh et al. (2012) PloS ONE 7: e44751.\n\
    iHS: BF Voight et al. (2006) PLoS Biology 4: e72.\n\
    EHH: PC Sabeti et al. (2002) Nature 419: 832–837.\n\
    \n\
    To normalize selscan output across frequency bins:\n\
    \n\
    selscan norm [--ihs|--xpehh|--nsl|--xpnsl|--ihh12] --files <file1.*.out> ... <fileN.*.out>\n\
    \n\
    To normalize selscan output and analyze non-overlapping windows of fixed bp for \n\
    extreme scores:\n\
    \n\
    selscan norm [--ihs|--xpehh|--nsl|--xpnsl|--ihh12] --files <file1.*.out> ... <fileN.*.out> --bp-win\n";

    const string ARG_FREQ_BINS = "--bins";
    const int DEFAULT_FREQ_BINS = 100;
    const string HELP_FREQ_BINS = "The number of frequency bins in [0,1] for score normalization.";

    const string ARG_FILES = "--files";
    const string DEFAULT_FILES = "__infile__";
    const string HELP_FILES = "A list of selscan output files, delimited by whitespace,\n\
    \tfor joint normalization. Each file starts with the header selscan writes:\n\
    \tiHS:    chr id pos freq ihh1 ihh0 ihs (four more columns with --ihs-detail)\n\
    \tnSL:    chr id pos freq ihh1 ihh0 nsl\n\
    \tXP-EHH: chr id pos gpos p1 ihh1 p2 ihh2 xpehh\n\
    \tXP-nSL: chr id pos gpos p1 sL1 p2 sL2 xpnsl\n\
    \tiHH12:  chr id pos freq ihh12\n\
    \tOutput from selscan before v3 has no header and no chr column, and is\n\
    \tnot accepted.";

    const string ARG_LOG = "--log";
    const string DEFAULT_LOG = "norm.log";
    const string HELP_LOG = "The log file name.";

    const string ARG_WINSIZE = "--winsize";
    const int DEFAULT_WINSIZE = 100000;
    const string HELP_WINSIZE = "The non-overlapping window size for calculating the percentage\n\
    \tof extreme SNPs.";


    const string ARG_WINOFFSET= "--win-offset";
    const int DEFAULT_WINOFFSET = 0;
    const string HELP_WINOFFSET = "The starting offset for non-overlapping windows. Defaults to 0.";

    const string ARG_QBINS = "--qbins";
    const int DEFAULT_QBINS = 10;
    const string HELP_QBINS = "Outlying windows are binned by number of sites within each\n\
    \twindow.  This is the number of quantile bins to use.";

    const string ARG_MINSNPS = "--min-snps";
    const int DEFAULT_MINSNPS = 10;
    const string HELP_MINSNPS = "Only consider a bp window if it has at least this many SNPs.";

    const string ARG_BPWIN = "--bp-win";
    const bool DEFAULT_BPWIN = false;
    const string HELP_BPWIN = "Use fixed-size bp windows (variable SNP count); outputs .windows with max/min score.\n\
    \tEach file must hold one chromosome, sorted by position.\n";

    const string ARG_IHS = "--ihs";
    const bool DEFAULT_IHS = false;
    const string HELP_IHS = "Do iHS normalization.";

    const string ARG_NSL = "--nsl";
    const bool DEFAULT_NSL = false;
    const string HELP_NSL = "Do nSL normalization.";

    const string ARG_XPEHH = "--xpehh";
    const bool DEFAULT_XPEHH = false;
    const string HELP_XPEHH = "Do XP-EHH normalization.";

    const string ARG_XPNSL = "--xpnsl";
    const bool DEFAULT_XPNSL = false;
    const string HELP_XPNSL = "Do XP-nSL normalization.";

    const string ARG_SOFT = "--ihh12";
    const bool DEFAULT_SOFT = false;
    const string HELP_SOFT = "Do ihh12 normalization.";

    const string ARG_FIRST = "--first";
    const bool DEFAULT_FIRST = false;
    const string HELP_FIRST = "Output only the first file's normalization.";

    const string ARG_CRIT_NUM = "--crit-val";
    const double DEFAULT_CRIT_NUM = 2;
    const string HELP_CRIT_NUM = "Set the critical value such that a SNP with |iHS| >= CRIT_VAL is marked as an extreme SNP.  Default as in Voight et al.";

    const string ARG_CRIT_PERCENT = "--crit-percent";
    const double DEFAULT_CRIT_PERCENT = -1;
    const string HELP_CRIT_PERCENT = "Set the critical value such that a SNP with iHS in the most extreme CRIT_PERCENT tails (two-tailed) is marked as an extreme SNP.\n\
    \tNot used by default.";

    // Added in v3
    const string ARG_LOG_INPUT = "--log-input";
    const string DEFAULT_LOG_INPUT = "__logfile__";
    const string HELP_LOG_INPUT =     "Specifies the log file used as an input for normalization.\n"
    "If provided, frequency-bin or mean/variance normalization is applied from the log-input file.\n"
    "Cannot be used together with --bins.\n"
    "Default: __logfile__.\n";

    // const string ARG_LOG_OUTPUT = "--log-output-only";
    // const string DEFAULT_LOG_OUTPUT = "__logfile__";
    // const string HELP_LOG_OUTPUT =     "Don't perform normalization, just output .\n"
    // "If provided, frequency-bin or mean/variance normalization is applied from the log-input file.\n"
    // "Cannot be used together with --bins.\n"
    // "Default: __logfile__.\n";


    const string ARG_BED = "--gene-bed"; // CHR START END GENE 
    const string DEFAULT_BED = "__filebed__"; 
    const string HELP_BED = "Provide a .bed file (<chr> <start> <end> <gene>) with gene annotations. ";


    const string ARG_GTF = "--gene-gtf";
    const string DEFAULT_GTF = "__filegtf__";
    const string HELP_GTF = "Provide a .gtf file with gene annotations. Extracts protein-coding genes only. Supports gzipped files";

    const string ARG_WIN_FILES = "--win-files";
    const string DEFAULT_WIN_FILES = "__filewin__";
    const string HELP_WIN_FILES = "Provide a .windows file generated by selscan norm to be annotated "
                                      "with gene names (requires --gene-bed).";

    const string ARG_NORM_FILES = "--norm-files";
    const string DEFAULT_NORM_FILES = "__filenorm__";
    const string HELP_NORM_FILES = "Provide a list of normalized selscan output files to be used for "
                                      "window-based output and/or gene annotation (requires --gene-bed and/or --bp-win).";

    const string ARG_GENE_OUT = "--gene-out";
    const string DEFAULT_GENE_OUT = "__geneout__";
    const string HELP_GENE_OUT = "File name for the gene table. By default it is written next to the first\n\
\tnormalized input, as <input>.<BED or GTF file name>.genetable.";

    const string ARG_GENE_SETA = "--gene-target";
    const string DEFAULT_GENE_SETA = "__genetarget__";
    const string HELP_GENE_SETA = "Provide a .genetable for target gene set for permutation tests.";

    const string ARG_GENE_SETB = "--gene-background";
    const string DEFAULT_GENE_SETB = "__genebackground__";
    const string HELP_GENE_SETB = "Provide a .genetable for background gene set for permutation test.";

    const string ARG_PERM_SEED = "--perm-seed";
    const int DEFAULT_PERM_SEED = -1;
    const string HELP_PERM_SEED = "Random seed for the gene-set permutation test (default: seeded from the clock; the seed used is printed).";

    const string ARG_FINE_PERCENTILE = "--fine-percentile";
    const bool DEFAULT_FINE_PERCENTILE = false;
    const string HELP_FINE_PERCENTILE = "If set, will use fine grain percentiles (1,2,3,...,100) for normalization.";

    const string ARG_LOG_ONLY = "--log-only";   
    const bool DEFAULT_LOG_ONLY = false;
    const string HELP_LOG_ONLY = "If set, outputs only the log file with normalization info (mean, variance) and skips normalized output. Useful with --log-input.";

    
    // const string ARG_NO_HEADER = "--no-header";
    // const bool DEFAULT_NO_HEADER = false;
    // const string HELP_NO_HEADER = "If set, will not include header in output files.";

    const int MISSING = -9999;

    // B18: numeric columns are read as tokens and parsed here, so a token
    // such as "nan", "inf" or "x" cannot put the stream into a fail state
    // (which silently zeroed every later row). readNum returns NaN for
    // anything that is not a finite number; readScore maps that to MISSING
    // and counts it in nBadRows (read pass only).
    long nBadRows = 0;
    static double parseNum(const string &t){
        char *e = nullptr;
        double v = strtod(t.c_str(), &e);
        if (t.empty() || e == t.c_str() || *e != '\0' || !std::isfinite(v)) return NAN;
        return v;
    }
    static double readNum(istream &in){ string t; in >> t; return parseNum(t); }
    double readScore(istream &in, bool count){
        double v = readNum(in);
        if (std::isnan(v)) { if (count) nBadRows++; return MISSING; }
        return v;
    }
    void reportBadRows(){
        if (nBadRows == 0) return;
        cerr << "WARNING: " << nBadRows << " rows have a non-numeric or non-finite score or frequency (e.g. nan, inf); they are treated as missing and not written.\n";
        flog << "WARNING: " << nBadRows << " rows have a non-numeric or non-finite score or frequency (e.g. nan, inf); they are treated as missing and not written.\n";
    }

    bool FINE_PERCENTILE = false;

    //returns number of lines in file
    //throws 0 if the file fails
    static vector<string> expectedHeaders(const string &stat, bool norm);
    int checkStatFile(ifstream &fin, const string &filename, const string &stat, bool norm);

    void readAllIHS(vector<string> filename, int fileLoci[], int nfiles, double freq[], double score[]);
    void readAllXPEHH(vector<string> filename, int fileLoci[], int nfiles, double freq1[], double freq2[], double score[]);
    void readAllIHH12(vector<string> filename, int fileLoci[], int nfiles, double freq1[], double score[]);

    // Returns the number of loci left out because their frequency (1.0) has no
    // bin; always 0 when includeFreqOne is true.
    int getMeanVarBins(double freq[], double data[], int nloci, double mean[], double variance[], int n[], int numBins, double threshold[], bool includeFreqOne);
    void windowQuantiles(const double sortedSNPs[], double frac[], int numWindows, int numQuantiles,
                         double quantileBound[], map<string, double> topWindowBoundary[], const string &label);
    double windowPercentile(int nSNPs, double frac, const double quantileBound[], int numQuantiles,
                            map<string, double> topWindowBoundary[]);
    int countLogBins(const std::string &binFile);
    void critCutoffs(const double score[], int nloci, double critPercent, double critNum,
                     bool oneTailed, double &upperCutoff, double &lowerCutoff);

    void normalizeIHSDataByBins(string &filename, string &outfilename, int &fileLoci, double mean[], double variance[], int n[], int numBins, double threshold[], double upperCutoff, double lowerCutoff, bool NSL);
    void normalizeXPEHHDataByBins(string &filename, string &outfilename, int &fileLoci, double mean[], double variance[], int n[], int numBins, double threshold[], double upperCutoff, double lowerCutoff, bool XPNSL);
    void normalizeIHH12DataByBins(string &filename, string &outfilename, int &fileLoci, double mean[], double variance[], int n[], int numBins, double threshold[], double upperCutoff, double lowerCutoff);

    void analyzeIHSBPWindows(string normedfiles[], int fileLoci[], int nfiles, int winSize, int numQuantiles, int minSNPs);
    void analyzeXPEHHBPWindows(string normedfiles[], int fileLoci[], int nfiles, int winSize, int numQuantiles, int minSNPs);
    void analyzeIHH12BPWindows(string normedfiles[], int fileLoci[], int nfiles, int winSize, int numQuantiles, int minSNPs);

    int countCols(ifstream &fin);
    int colsToSkip(ifstream &fin, int numCols);
    void skipCols(ifstream &fin, int numCols);

    int countFields(const string &str);
    bool isint(string str);

    //------------

    ofstream flog;
    int runToolNorm(int argc, char *argv[]);

    void getMeanVarBinsFromLog(const std::string &binFile,
                                       double freq[], double data[], int nloci,
                                       double mean[], double variance[], int n[],
                                       int numBins, double threshold[], bool XPORSOFT);
    
        
    double mean(const std::vector<double> &v) {
        return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
    }

};

// int main(int argc, char *argv[]){
//     std::vector<double> lengths = {1000, 2000, 1500, 3000};
//     std::vector<double> scores = {2.5, 3.0, 2.8, 3.5};
//     std::vector<double> residuals = SelscanNorm::regress_out_length(lengths, scores);
//     for (double res : residuals) {
//         std::cout << res << std::endl;
//     }

//     return SelscanNorm::runToolNorm(argc, argv);
// }

#endif


