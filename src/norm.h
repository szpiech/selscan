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

#include <iostream>
#include <fstream>
#include <map>
#include <cmath>
#include <cstdlib>
#include <cctype>
#include <vector>
#include <cstring>
#include <sstream>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include "param_t.h"

using namespace std;


// added afterwards
#include <gsl/gsl_multifit.h>

struct BinData {
    std::vector<double> bins;
    std::vector<double> means;
    std::vector<double> variances;
};




class SelscanNorm{
    public:
    const string VERSION = "1.3.1";

    const string PREAMBLE = " -- a program for downstream analysis of selscan output\n\
    Source code and binaries can be found at\n\
    \t<https://www.github.com/szpiech/selscan>\n\
    \n\
    Citations:\n\
    \n\
    selscan: ZA Szpiech and RD Hernandez (2014) MBE, 31: 2824-2827.\n\
    iHH12: R Torres et al. (2018) PLoS Genetics 15: e1007898.\n\
        N Garud et al. (2015) PLoS Genetics, 11: 1–32.\n\
    nSL: A Ferrer-Admetlla, et al. (2014) MBE, 31: 1275-1291.\n\
    XP-nSL: Szpiech et al. (2020) bioRxiv doi: \n\
            https://doi.org/10.1101/2020.05.19.104380.\n\
    XP-EHH: PC Sabeti et al. (2007) Nature, 449: 913–918.\n\
    iHS: BF Voight et al. (2006) PLoS Biology, 4: e72.\n\
    \n\
    To normalize selscan output across frequency bins:\n\
    \n\
    ./norm [--ihs|--xpehh|--nsl|--xpnsl|--ihh12] --files <file1.*.out> ... <fileN.*.out>\n\
    \n\
    To normalize selscan output and analyze non-overlapping windows of fixed bp for \n\
    extreme scores:\n\
    \n\
    ./norm [--ihs|--xpehh|--nsl|--xpnsl|--ihh12] --files <file1.*.out> ... <fileN.*.out> --bp-win\n";

    const string ARG_FREQ_BINS = "--bins";
    const int DEFAULT_FREQ_BINS = 100;
    const string HELP_FREQ_BINS = "The number of frequency bins in [0,1] for score normalization.";

    const string ARG_FILES = "--files";
    const string DEFAULT_FILES = "infile";
    const string HELP_FILES = "A list of files delimited by whitespace for\n\
    \tjoint normalization.\n\
    \tExpected format for iHS or nSL files (no header):\n\
    \t<locus name> <physical pos> <freq> <ihh1/sL1> <ihh2/sL2> <ihs/nsl>\n\
    \tExpected format for XP-EHH files (one line header):\n\
    \t<locus name> <physical pos> <genetic pos> <freq1> <ihh1> <freq2> <ihh2> <xpehh>\n\
    \tExpected format for iHH12 files (one line header):\n\
    \t<locus name> <physical pos> <freq1> <ihh12>";

    const string ARG_LOG = "--log";
    const string DEFAULT_LOG = "logfile";
    const string HELP_LOG = "The log file name.";

    const string ARG_WINSIZE = "--winsize";
    const int DEFAULT_WINSIZE = 100000;
    const string HELP_WINSIZE = "The non-overlapping window size for calculating the percentage\n\
    \tof extreme SNPs.";

    const string ARG_QBINS = "--qbins";
    const int DEFAULT_QBINS = 10;
    const string HELP_QBINS = "Outlying windows are binned by number of sites within each\n\
    \twindow.  This is the number of quantile bins to use.";

    const string ARG_MINSNPS = "--min-snps";
    const int DEFAULT_MINSNPS = 10;
    const string HELP_MINSNPS = "Only consider a bp window if it has at least this many SNPs.";

    const string ARG_SNPWIN = "--snp-win";
    const bool DEFAULT_SNPWIN = false;
    const string HELP_SNPWIN = "<not implemented> If set, will use windows of a constant\n\
    \tSNP size with varying bp length.";

    const string ARG_SNPWINSIZE = "--snp-win-size";
    const int DEFAULT_SNPWINSIZE = 50;
    const string HELP_SNPWINSIZE = "<not implemented> The number of SNPs in a window.";

    const string ARG_BPWIN = "--bp-win";
    const bool DEFAULT_BPWIN = false;
    const string HELP_BPWIN = "If set, will use windows of a constant bp size with varying\n\
    \tnumber of SNPs.";

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
    const string HELP_CRIT_NUM = "Set the critical value such that a SNP with |iHS| > CRIT_VAL is marked as an extreme SNP.  Default as in Voight et al.";

    const string ARG_CRIT_PERCENT = "--crit-percent";
    const double DEFAULT_CRIT_PERCENT = -1;
    const string HELP_CRIT_PERCENT = "Set the critical value such that a SNP with iHS in the most extreme CRIT_PERCENT tails (two-tailed) is marked as an extreme SNP.\n\
    \tNot used by default.";

    // Added in v3
    const string ARG_LOG_INPUT = "--log-input";
    const string DEFAULT_LOG_INPUT = "logfile";
    const string HELP_LOG_INPUT = "The log file name.";

    const string ARG_BED = "--gene-bed";
    const string DEFAULT_BED = "__bedfile__";
    const string HELP_BED = "Provide a BED file to get annotations for outlier windows.";

    const string ARG_FINE_PERCENTILE = "--fine-perc";
    const bool DEFAULT_FINE_PERCENTILE = false;
    const string HELP_FINE_PERCENTILE = "If set, will use fine grain percentiles (1,2,3,...,100) for normalization.";

    const int MISSING = -9999;
    bool FINE_PERCENTILE = false;

        //returns number of lines in file
        //throws 0 if the file fails
        int checkIHSfile(ifstream &fin);
        int checkXPEHHfile(ifstream &fin);
        int checkIHH12file(ifstream &fin);

        void readAllIHS(vector<string> filename, int fileLoci[], int nfiles, double freq[], double score[]);
        void readAllXPEHH(vector<string> filename, int fileLoci[], int nfiles, double freq1[], double freq2[], double score[]);
        void readAllIHH12(vector<string> filename, int fileLoci[], int nfiles, double freq1[], double score[]);

        void getMeanVarBins(double freq[], double data[], int nloci, double mean[], double variance[], int n[], int numBins, double threshold[]);

        void normalizeIHSDataByBins(string &filename, string &outfilename, int &fileLoci, double mean[], double variance[], int n[], int numBins, double threshold[], double upperCutoff, double lowerCutoff);
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

    ofstream flog;
    int runToolNorm(int argc, char *argv[]);
    std::vector<double> regress_out_length(const std::vector<double> &lengths,
                                       const std::vector<double> &scores);



    void getMeanVarBinsFromLog(const std::string &binFile,
                                       double freq[], double data[], int nloci,
                                       double mean[], double variance[], int n[],
                                       int numBins, double threshold[]);
    // BinData getMeanVarsFromLog(const std::string &filename) {
    //     std::ifstream infile(filename);
    //     if (!infile.is_open()) {
    //         throw std::runtime_error("Error opening file: " + filename);
    //     }

    //     BinData data;
    //     std::string line;

    //     while (std::getline(infile, line)) {
    //         if (line.empty()) continue;
    //         if (line.rfind("Total", 0) == 0) continue;  // skip "Total loci..."
    //         if (line.rfind("bin", 0) == 0) continue;    // skip header

    //         std::istringstream iss(line);
    //         double bin, mean, var;
    //         int num;

    //         iss >> bin >> num >> mean >> var;

    //         // Skip if num == 0 or mean/var are NaN placeholders
    //         if (num == 0 || !iss) continue;

    //         data.bins.push_back(bin);
    //         data.means.push_back(mean);
    //         data.variances.push_back(var);
    //     }

    //     return data;
    // }                                 
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