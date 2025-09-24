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

#include <gsl/gsl_multifit.h>

using namespace std;

const string VERSION = "2.0.0";

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

///////-------------
const string ARG_LOG_INPUT = "--log-input";
const string DEFAULT_LOG_INPUT = "logfile";
const string HELP_LOG_INPUT = "The log file name.";


const string ARG_BED = "--gene-bed";
const string DEFAULT_BED = "__bedfile__";
const string HELP_BED = "Provide a BED file to get annotations for outlier windows.";

// const string ARG_FIXED_SNPS = "--fixed-snps";
// const int DEFAULT_FIXED_SNPS = 100;
// const std::string HELP_FIXED_SNPS = "Fix the number of SNPs per window/bin instead of\n\
// \tusing a fixed window size or quantile bins. The argument specifies the\n\
// \tnumber of SNPs per bin.";
// each bin to represent roughly the same number of SNPs, not the same number of windows.

//Equal number of windows per bin (by rank)
//Each bin has approx fixed SNP count total
//////--------------

const string ARG_WINSIZE = "--winsize";
const int DEFAULT_WINSIZE = 100000;
const string HELP_WINSIZE = "The non-overlapping window size for calculating the percentage\n\
\tof extreme SNPs.";
//Windows are grouped by rank, so bins always have equal numbers of windows, no matter how SNP counts vary.

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

const int MISSING = -9999;


enum StatType { IHS, XPEHH, IHH12 };

ofstream flog;
vector<int> fileLoci;
vector<string> outfilename;
int totalLoci = 0;


// ---- Configuration extracted from params ----
struct Config {
    int numBins;
    vector<string> filenames;
    int winSize;
    string infoOutfile;
    int numQBins;
    int minSNPs;
    int snpWinSize;
    bool BPWIN;
    bool SNPWIN;
    bool FIRST;
    double critNum;
    double critPercent;
    bool IHS;
    bool NSL;
    bool SOFT;
    bool XPEHH;
    bool XPNSL;
};


void getMeanVarBins(double freq[], double data[], int nloci, double mean[], double variance[], int n[], int numBins, double threshold[], string normLogInput = "");
void normalizeDataByBins(const string& filename,
                     const string& outfilename,
                     const int& fileLoci,
                     double mean[], double variance[], int n[],
                     int numBins, double threshold[],
                     double upperCutoff, double lowerCutoff,
                     const string& type,
                     bool XPNSL = false);

void analyzeBPWindows(const string normedfiles[], const int fileLoci[], int nfiles, int winSize, int numQuantiles, int minSNPs, StatType statType);


int checkFile( std::ifstream &fin, const std::string &name);
void readAll(const vector<string>& filenames,
             const int fileLoci[],
             int nfiles,
             const string& type,
             double* freq1,
             double* score,
             double* freq2 = nullptr);

static void chooseCutoffs(const Config& c, double* scoreSorted, int totalLoci,
                          double& upperCutoff, double& lowerCutoff);

static bool validateConfig(const Config& c) {
    if (c.numBins <= 0) {
        cerr << "ERROR: Must have a positive integer of frequency bins.\n";
        return false;
    }
    if (c.numQBins <= 0) {
        cerr << "ERROR: Must have a positive integer of quantile bins.\n";
        return false;
    }
    if (c.winSize <= 0) {
        cerr << "ERROR: Must have a positive integer window size.\n";
        return false;
    }
    if (c.critNum <= 0) {
        cerr << "ERROR: Must give a positive critical value for |iHS| scores.\n";
        return false;
    }
    if (c.critPercent != DEFAULT_CRIT_PERCENT && (c.critPercent <= 0 || c.critPercent >= 1)) {
        cerr << "ERROR: Critical percentile must be in (0,1).\n";
        return false;
    }
    int chosen = (c.IHS?1:0) + (c.XPEHH?1:0) + (c.NSL?1:0) + (c.SOFT?1:0);
    if (chosen != 1) {
        cerr << "Must specify exactly one of " << ARG_IHS << ", " << ARG_XPEHH
             << ", " << ARG_NSL << ", " << ARG_SOFT << ", " << ARG_XPNSL << ".\n";
        return false;
    }
    return true;
}


static bool openLog(const Config& c, int argc, char* argv[]) {
    flog.open(c.infoOutfile.c_str());
    if (flog.fail()) {
        cerr << "ERROR: " << c.infoOutfile << " " << strerror(errno) << endl;
        return false;
    }
    for (int i = 0; i < argc; i++) flog << argv[i] << " ";
    flog << "\n\n";
    return true;
}


static bool prepareFiles(const Config& c,
                         vector<string>& outfilename,
                         vector<int>& fileLoci,
                         int& totalLoci) {

                            //fileLoci: number of loci per input file
    const int nfiles = static_cast<int>(c.filenames.size());
    outfilename.resize(nfiles);
    fileLoci.resize(nfiles);
    totalLoci = 0;

    for (int i = 0; i < nfiles; i++) {
        char str[16];
        sprintf(str, "%d", c.numBins);
        if (c.IHS || c.NSL) outfilename[i] = c.filenames[i] + "." + string(str) + "bins.norm";
        if (c.XPEHH || c.SOFT) outfilename[i] = c.filenames[i] + ".norm";

        ifstream fin(c.filenames[i].c_str());
        if (fin.fail()) {
            cerr << "ERROR: " << c.infoOutfile << " " << strerror(errno);
            flog << "ERROR: " << c.infoOutfile << " " << strerror(errno);
            return false;
        } else {
            cerr << "Opened " << c.filenames[i] << endl;
        }

        try {
            if (c.IHS || c.NSL) fileLoci[i] = checkFile(fin, "ihs");
            if (c.SOFT)         fileLoci[i] = checkFile(fin, "ihh12");
            if (c.XPEHH)        fileLoci[i] = checkFile(fin, "xpehh");
            totalLoci += fileLoci[i];
        } catch (...) {
            return false;
        }
        fin.close();
    }

    cerr << "\nTotal loci: " << totalLoci << endl;
    flog << "\nTotal loci: " << totalLoci << endl;
    return true;
}


static void computeThresholdsAndPrintBins(int numBins,
                                          double* mean, double* variance, int* n,
                                          const double* threshold,
                                          bool printFreqThresholds) {
    if (printFreqThresholds) {
        cerr << "bin\tnum\tmean\tvariance\n";
        flog << "bin\tnum\tmean\tvariance\n";
        for (int i = 0; i < numBins; i++) {
            cerr << threshold[i] << "\t" << n[i] << "\t" << mean[i] << "\t" << variance[i] << endl;
            flog << threshold[i] << "\t" << n[i] << "\t" << mean[i] << "\t" << variance[i] << endl;
        }
    } else {
        cerr << "num\tmean\tvariance\n";
        flog << "num\tmean\tvariance\n";
        for (int i = 0; i < numBins; i++) {
            cerr << n[i] << "\t" << mean[i] << "\t" << variance[i] << endl;
            flog << n[i] << "\t" << mean[i] << "\t" << variance[i] << endl;
        }
    }
}




static void processIHSorNSL(const Config& c,
                            const vector<string>& outfilename,
                            const vector<int>& fileLoci,
                            int totalLoci) {
    cerr << "Reading all data.\n";
    double* freq  = new double[totalLoci];
    double* score = new double[totalLoci];

    readAll(c.filenames, fileLoci.data(), (int)c.filenames.size(), "ihs", freq, score);

    // Prepare bins (uniform over [0,1])
    double minFreq = 0.0, maxFreq = 1.0;
    double step = (maxFreq - minFreq) / double(c.numBins);
    double* threshold = new double[c.numBins];
    for (int b = 0; b < c.numBins; b++) threshold[b] = minFreq + (b + 1) * step;

    double* mean     = new double[c.numBins];
    double* variance = new double[c.numBins];
    int*    n        = new int[c.numBins];

    cerr << "Calculating mean and variance per frequency bin:\n\n";
    getMeanVarBins(freq, score, totalLoci, mean, variance, n, c.numBins, threshold);

    // cutoffs
    gsl_sort(score, 1, totalLoci);
    double upperCutoff, lowerCutoff;
    chooseCutoffs(c, score, totalLoci, upperCutoff, lowerCutoff);

    delete [] freq;
    delete [] score;

    // Print bins info
    computeThresholdsAndPrintBins(c.numBins, mean, variance, n, threshold, /*printFreqThresholds=*/true);

    // Normalize per file
    int nfiles = (c.FIRST ? 1 : (int)c.filenames.size());
    for (int i = 0; i < nfiles; i++) {
        cerr << "Normalizing " << c.filenames[i] << "\n";
        normalizeDataByBins(c.filenames[i], outfilename[i], fileLoci[i],
                               mean, variance, n, c.numBins, threshold,
                               upperCutoff, lowerCutoff, "ihs");
    }

    // cleanup
    delete [] threshold;
    delete [] mean;
    delete [] variance;
    delete [] n;

    // Window analyses
    if (c.BPWIN) {
        analyzeBPWindows(const_cast<string*>(outfilename.data()), fileLoci.data(),
                            nfiles, c.winSize, c.numQBins, c.minSNPs, IHS);
    }
    // if (c.SNPWIN) analyzeSNPWindows(...); // original commented-out line
}


static void processXPEHHorSOFT(const Config& c,
                               const vector<string>& outfilename,
                               const vector<int>& fileLoci,
                               int totalLoci) {
    cerr << "Reading all data.\n";
    double* freq1 = new double[totalLoci];
    double* freq2 = nullptr;
    if (c.XPEHH) freq2 = new double[totalLoci];
    double* score = new double[totalLoci];

    if (c.XPEHH) readAll(c.filenames, fileLoci.data(), (int)c.filenames.size(), "xpehh", freq1, score, freq2);
    if (c.SOFT)  readAll(c.filenames, fileLoci.data(), (int)c.filenames.size(), "ihh12", freq1, score);

    // For these modes, force a single bin:
    int numBins = 1;
    double minFreq = 0.0, maxFreq = 1.0;
    double step = (maxFreq - minFreq) / double(numBins);
    double* threshold = new double[numBins];
    for (int b = 0; b < numBins; b++) threshold[b] = minFreq + (b + 1) * step;

    double* mean     = new double[numBins];
    double* variance = new double[numBins];
    int*    n        = new int[numBins];

    cerr << "Calculating mean and variance:\n\n";
    getMeanVarBins(freq1, score, totalLoci, mean, variance, n, numBins, threshold);

    gsl_sort(score, 1, totalLoci);
    double upperCutoff, lowerCutoff;
    chooseCutoffs(c, score, totalLoci, upperCutoff, lowerCutoff);

    delete [] freq1;
    if (c.XPEHH) delete [] freq2;
    delete [] score;

    // Print bins info
    computeThresholdsAndPrintBins(numBins, mean, variance, n, threshold, /*printFreqThresholds=*/false);

    // Normalize per file
    int nfiles = (c.FIRST ? 1 : (int)c.filenames.size());
    for (int i = 0; i < nfiles; i++) {
        cerr << "Normalizing " << c.filenames[i] << "\n";
        if (c.XPEHH) {
            normalizeDataByBins(c.filenames[i], outfilename[i], fileLoci[i],
                                     mean, variance, n, numBins, threshold,
                                     upperCutoff, lowerCutoff, "xpehh", c.XPNSL);
        }
        if (c.SOFT) {
            normalizeDataByBins(c.filenames[i], outfilename[i], fileLoci[i],
                                     mean, variance, n, numBins, threshold,
                                     upperCutoff, lowerCutoff, "ihh12");
        }
    }

    // cleanup
    delete [] threshold;
    delete [] mean;
    delete [] variance;
    delete [] n;

    // Window analyses
    if (c.BPWIN) {
        if (c.XPEHH) analyzeBPWindows(const_cast<string*>(outfilename.data()), fileLoci.data(),
                                           nfiles, c.winSize, c.numQBins, c.minSNPs, XPEHH);
        if (c.SOFT)  analyzeBPWindows(const_cast<string*>(outfilename.data()), fileLoci.data(),
                                           nfiles, c.winSize, c.numQBins, c.minSNPs, IHH12);
    }else if (c.SNPWIN) {
        cerr << "\nSNP window analysis not yet implemented for XPEHH or IHH12.\n\n";
    }
}


struct Stats {
    long long n;   // sample size
    double mean;   // mean
    double var;    // variance
};

// Combine stats across files for a single bin
Stats combine_stats(const vector<Stats>& all) {
    long long total_n = 0;
    double weighted_sum = 0.0;

    for (auto& s : all) {
        total_n += s.n;
        weighted_sum += s.n * s.mean;
    }
    double global_mean = weighted_sum / total_n;

    double numerator = 0.0;
    for (auto& s : all) {
        numerator += (s.n - 1) * s.var;                 // within-file
        numerator += s.n * (s.mean - global_mean) * (s.mean - global_mean); // between-file
    }
    double global_var = numerator / (total_n - 1);

    return { total_n, global_mean, global_var };
}

// Combine stats per bin across multiple files
vector<Stats> combine_per_bin(const vector<vector<Stats>>& files_per_bin) {
    if (files_per_bin.empty()) return {};
    size_t num_bins = files_per_bin[0].size();

    vector<Stats> result(num_bins);

    for (size_t b = 0; b < num_bins; b++) {
        vector<Stats> bin_stats;
        for (auto& file : files_per_bin) {
            bin_stats.push_back(file[b]);
        }
        result[b] = combine_stats(bin_stats);
    }
    return result;
}


struct BedEntry {
    int start;
    int end;
    std::string gene; // gene name or ID
};


class Bin {
public:
    int num;        // number of observations
    double mean;
    double variance;

    // Constructor for empty bin
    Bin() : num(0), mean(0.0), variance(0.0) {}

    // Method to update bin with a new observation
    void update(double x) {
        num += 1;

        if (num == 1) {
            // first observation
            mean = x;
            variance = 0.0;
        } else {
            double delta = x - mean;
            mean += delta / num;
            double delta2 = x - mean;
            variance = ((num - 2) * variance + delta * delta2) / (num - 1);
        }
    }
};



std::vector<BedEntry> readAndSortBedForChrom(const std::string &filename, const std::string &targetChrom) {
    std::ifstream in(filename);
    if (!in) {
        throw std::runtime_error("Cannot open BED file: " + filename);
    }

    std::vector<BedEntry> entries;
    std::string line, chrom;
    int start, end;

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        iss >> chrom >> start >> end;

        if (chrom != targetChrom) continue; // skip unwanted chromosomes

        std::string extra;
        std::getline(iss, extra); // remaining columns
        entries.push_back({start, end, extra});
    }

    std::sort(entries.begin(), entries.end(), [](const BedEntry &a, const BedEntry &b) {
        return a.start < b.start;
    });

    return entries;
}



static void chooseCutoffs(const Config& c, double* scoreSorted, int totalLoci,
                          double& upperCutoff, double& lowerCutoff) {
    if (c.critPercent != DEFAULT_CRIT_PERCENT && (c.critPercent > 0 && c.critPercent < 1)) {
        upperCutoff = gsl_stats_quantile_from_sorted_data(scoreSorted, 1, totalLoci, 1 - c.critPercent / 2.0);
        lowerCutoff = gsl_stats_quantile_from_sorted_data(scoreSorted, 1, totalLoci, c.critPercent / 2.0);

        cerr << "\nTop cutoff: " << upperCutoff << endl;
        cerr << "Bottom cutoff: " << lowerCutoff << "\n\n";
        flog << "\nTop cutoff: " << upperCutoff << endl;
        flog << "Bottom cutoff: " << lowerCutoff << "\n\n";
    } else {
        upperCutoff = c.critNum;
        lowerCutoff = -c.critNum;
    }
}


int countCols(ifstream &fin);
int colsToSkip(ifstream &fin, int numCols);
void skipCols(ifstream &fin, int numCols);

int countFields(const string &str);
bool isint(string str);



// ---------- Helpers ----------
static bool parseParams(int argc, char* argv[], param_t& params) {
    try {
        params.setPreamble(PREAMBLE);
        params.addFlag(ARG_FREQ_BINS, DEFAULT_FREQ_BINS, "", "");
        params.addListFlag(ARG_FILES, DEFAULT_FILES, "", "");
        params.addFlag(ARG_LOG, DEFAULT_LOG, "", "");
        params.addFlag(ARG_WINSIZE, DEFAULT_WINSIZE, "", "");
        params.addFlag(ARG_QBINS, DEFAULT_QBINS, "", "");
        params.addFlag(ARG_MINSNPS, DEFAULT_MINSNPS, "", "");
        params.addFlag(ARG_SNPWIN, DEFAULT_SNPWIN, "SILENT", "");
        params.addFlag(ARG_SNPWINSIZE, DEFAULT_SNPWINSIZE, "SILENT", "");
        params.addFlag(ARG_BPWIN, DEFAULT_BPWIN, "", "");
        params.addFlag(ARG_FIRST, DEFAULT_FIRST, "", "");
        params.addFlag(ARG_CRIT_NUM, DEFAULT_CRIT_NUM, "", "");
        params.addFlag(ARG_CRIT_PERCENT, DEFAULT_CRIT_PERCENT, "", "");
        params.addFlag(ARG_IHS, DEFAULT_IHS, "", "");
        params.addFlag(ARG_NSL, DEFAULT_NSL, "", "");
        params.addFlag(ARG_SOFT, DEFAULT_SOFT, "", "");
        params.addFlag(ARG_XPEHH, DEFAULT_XPEHH, "", "");
        params.addFlag(ARG_XPNSL, DEFAULT_XPNSL, "", "");
        params.parseCommandLine(argc, argv);
        return true;
    } catch (...) {
        return false;
    }
}

bool fillConfigFromParams(param_t& p, Config& c) {
    c.numBins     = p.getIntFlag(ARG_FREQ_BINS);
    c.filenames   = p.getStringListFlag(ARG_FILES);
    c.winSize     = p.getIntFlag(ARG_WINSIZE);
    c.infoOutfile = p.getStringFlag(ARG_LOG);
    c.numQBins    = p.getIntFlag(ARG_QBINS);
    c.minSNPs     = p.getIntFlag(ARG_MINSNPS);
    c.snpWinSize  = p.getIntFlag(ARG_SNPWINSIZE);
    c.BPWIN       = p.getBoolFlag(ARG_BPWIN);
    c.SNPWIN      = p.getBoolFlag(ARG_SNPWIN);
    c.FIRST       = p.getBoolFlag(ARG_FIRST);
    c.critNum     = p.getDoubleFlag(ARG_CRIT_NUM);
    c.critPercent = p.getDoubleFlag(ARG_CRIT_PERCENT);
    c.IHS         = p.getBoolFlag(ARG_IHS);
    c.NSL         = p.getBoolFlag(ARG_NSL);
    c.SOFT        = p.getBoolFlag(ARG_SOFT);
    c.XPEHH       = p.getBoolFlag(ARG_XPEHH);
    c.XPNSL       = p.getBoolFlag(ARG_XPNSL);

    if (c.XPNSL) c.XPEHH = true; // original behavior
    return true;
}

static inline double quantile_from_sorted(const std::vector<double>& v, double p) {
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    if (p <= 0) return v.front();
    if (p >= 1) return v.back();
    double idx = p * (v.size() - 1);
    size_t lo = (size_t)std::floor(idx);
    size_t hi = std::min(lo + 1, v.size() - 1);
    double w = idx - lo;
    return (1.0 - w) * v[lo] + w * v[hi];
}

struct SNP_IHS {
    int    pos;
    double normed;
    int    crit;      // 0/1
};
struct SNP_XPEHH {
    int    pos;
    double normed;
    int    crit;      // -1/0/1
};
struct SNP_IHH12 {
    int    pos;
    int    crit;      // 0/1
};


// --------------------- Function: Validate parameters ---------------------
bool validateParams(param_t &params, int &numBins, int &numQBins, int &winSize,
                    double &critNum, double &critPercent, bool &IHS, bool &NSL,
                    bool &SOFT, bool &XPEHH, bool &XPNSL) 
{
    numBins = params.getIntFlag(ARG_FREQ_BINS);
    numQBins = params.getIntFlag(ARG_QBINS);
    winSize = params.getIntFlag(ARG_WINSIZE);
    critNum = params.getDoubleFlag(ARG_CRIT_NUM);
    critPercent = params.getDoubleFlag(ARG_CRIT_PERCENT);
    IHS = params.getBoolFlag(ARG_IHS);
    NSL = params.getBoolFlag(ARG_NSL);
    SOFT = params.getBoolFlag(ARG_SOFT);
    XPEHH = params.getBoolFlag(ARG_XPEHH);
    XPNSL = params.getBoolFlag(ARG_XPNSL);

    if(XPNSL) XPEHH = true;

    if(numBins <= 0) { std::cerr << "ERROR: Must have a positive integer of frequency bins.\n"; return false; }
    if(numQBins <= 0) { std::cerr << "ERROR: Must have a positive integer of quantile bins.\n"; return false; }
    if(winSize <= 0) { std::cerr << "ERROR: Must have a positive integer window size.\n"; return false; }
    if(critNum <= 0) { std::cerr << "ERROR: Must give a positive critical value for |iHS| scores.\n"; return false; }
    if(critPercent != DEFAULT_CRIT_PERCENT && (critPercent <= 0 || critPercent >= 1)) { 
        std::cerr << "ERROR: Critical percentile must be in (0,1).\n"; 
        return false; 
    }
    if(IHS + XPEHH + NSL + SOFT != 1){
        std::cerr << "Must specify exactly one of " + ARG_IHS + ", " + ARG_XPEHH + "," + ARG_NSL + "," + ARG_SOFT + "," + ARG_XPNSL + ".\n";
        return false;
    }
    return true;
}



// --------------------- Function: Log command line ---------------------
void logCommandLine(int argc, char *argv[], const std::string& infoOutfile) {
    flog.open(infoOutfile.c_str());
    if(flog.fail()){
        std::cerr << "ERROR: Cannot open log file " << infoOutfile << ": " << strerror(errno) << std::endl;
        exit(1);
    }
    for(int i = 0; i < argc; i++) flog << argv[i] << " ";
    flog << "\n\n";
}


// ---------- main ----------
int main(int argc, char* argv[]) {
    cerr << "norm v" << VERSION << "\n";

    // 1) Parse params
    param_t params;
    if (!parseParams(argc, argv, params)) return 1;

    // 2) Fill config + basic validation
    Config cfg;
    fillConfigFromParams(params, cfg);
    if (!validateConfig(cfg)) return 1;

    // 3) Logging
    if (!openLog(cfg, argc, argv)) return 1; // affect flog

    // 4) Prepare files and count loci
    cerr << "You have provided " << cfg.filenames.size()
         << " output files for joint normalization.\n";

    if (!prepareFiles(cfg, outfilename, fileLoci, totalLoci)) return 1;

    // 5) Branch: IHS/NSL vs XPEHH/SOFT
    if (cfg.IHS || cfg.NSL) {
        processIHSorNSL(cfg, outfilename, fileLoci, totalLoci);
    } else if (cfg.XPEHH || cfg.SOFT) {
        processXPEHHorSOFT(cfg, outfilename, fileLoci, totalLoci);
    }

    // 6) Done
    flog.close();
    return 0;
}


// if any part of gene overlaps a window, mark that window as overlapping the gene
// if a window overlaps multiple genes, list all genes in the annotation (comma-separated?), 
/**
 * Regresses score ~ intercept + log(length) using GSL
 * and returns residuals (length-corrected scores).
 *
 * @param lengths  Vector of gene lengths
 * @param scores   Vector of per-gene scores
 * @return residuals vector (score corrected for length)
 */
std::vector<double> regress_out_length(const std::vector<double> &lengths,
                                       const std::vector<double> &scores) {
    size_t n = lengths.size();
    size_t p = 2; // intercept + log(length)

    gsl_matrix *X = gsl_matrix_alloc(n, p);
    gsl_vector *y = gsl_vector_alloc(n);
    gsl_vector *c = gsl_vector_alloc(p);
    gsl_matrix *cov = gsl_matrix_alloc(p, p);
    double chisq;

    // Fill design matrix and response
    for (size_t i = 0; i < n; i++) {
        gsl_matrix_set(X, i, 0, 1.0);                  // intercept
        gsl_matrix_set(X, i, 1, std::log(lengths[i])); // log(length)
        gsl_vector_set(y, i, scores[i]);
    }

    // Fit linear regression
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n, p);
    gsl_multifit_linear(X, y, c, cov, &chisq, work);
    gsl_multifit_linear_free(work);

    double beta0 = gsl_vector_get(c, 0);
    double beta1 = gsl_vector_get(c, 1);

    // Compute residuals
    std::vector<double> residuals(n);
    for (size_t i = 0; i < n; i++) {
        double pred = beta0 + beta1 * std::log(lengths[i]);
        residuals[i] = scores[i] - pred;
    }

    // Free memory
    gsl_matrix_free(X);
    gsl_matrix_free(cov);
    gsl_vector_free(y);
    gsl_vector_free(c);

    return residuals;
}



void analyzeBPWindows(const string normedfiles[], const int fileLoci[], int nfiles, int winSize, int numQuantiles, int minSNPs, StatType statType)
{
    cerr << "\nAnalyzing BP windows:\n\n";

    vector<int> *winStarts = new vector<int>[nfiles]; // Vector of vectors to hold the start positions of each window
    vector<int> *nSNPs = new vector<int>[nfiles]; // Vector of vectors to hold the number of SNPs in each window

    // For iHS and IHH12
    vector<double> *fracCrit = nullptr; // Fraction of critical SNPs in each window
    vector<double> *maxAbsScore = nullptr; // Maximum absolute score in each window

    // For XPEHH
    vector<double> *fracCritTop = nullptr; // Fraction of critical SNPs in the top tail
    vector<double> *fracCritBot = nullptr; // Fraction of critical SNPs in the bottom tail
    vector<double> *maxScore = nullptr; // Maximum score in each window
    vector<double> *minScore = nullptr; // Minimum score in each window

    if (statType == IHS || statType == IHH12)
    {
        fracCrit = new vector<double>[nfiles];
        maxAbsScore = new vector<double>[nfiles];
    }
    else if (statType == XPEHH)
    {
        fracCritTop = new vector<double>[nfiles];
        fracCritBot = new vector<double>[nfiles];
        maxScore = new vector<double>[nfiles];
        minScore = new vector<double>[nfiles];
    }

    ifstream fin;
    ofstream fout;
    string *winfilename = new string[nfiles]; // Output file names for windows

    // Create a string to hold the window size in kilobases
    // This is used to create the output file names
    // e.g., "100kb.windows" for a winSize of 100000
    // The size is divided by 1000 to convert from bp to kb
    // and stored as a string for the file name
    // This is done to avoid using floating-point numbers in file names
    // which can cause issues on some filesystems
    // The size is stored as an integer in kilobases
    // and formatted as a string with sprintf
    // to ensure it is always 3 digits long (e.g., "001", "010", "100")
    // This allows for consistent file naming and easy identification of window sizes


    char str[10];
    sprintf(str, "%d", winSize / 1000);

    // Variables used in all types (some unused in some)
    double gpos;
    string name, header;
    int pos;
    double freq, ihh1, ihh2, data, normedData;
    bool critBool;
    int critInt;

    // for xpehh
    double freq1, freq2;

    // for ihh12
    double ihh12;

    int numWindows = 0, numWindowsTop = 0, numWindowsBot = 0;

    for (int i = 0; i < nfiles; i++)
    {
        fin.open(normedfiles[i].c_str());
        if (fin.fail())
        {
            cerr << "ERROR: " << normedfiles[i] << " " << strerror(errno);
            throw -1;
        }

        if (statType == XPEHH || statType == IHH12)
        {
            getline(fin, header);  // Skip header line
        }

        //generate winfile names
        winfilename[i] = normedfiles[i];
        winfilename[i] += ".";
        winfilename[i] += str;
        winfilename[i] += "kb.windows";

        //Load information into vectors for analysis
        int winStart = 1;
        int winEnd = winStart + winSize - 1;
        int numSNPsWin = 0;

        // Initialize counters and extrema
        int numCrit = 0;
        int numCritTop = 0;
        int numCritBot = 0;
        double maxAbs = -99999.9;
        double maxVal = -99999.9;
        double minVal = 99999.9;

        for (int j = 0; j < fileLoci[i]; j++)
        {
            if (statType == IHS)
            {
                fin >> name >> pos >> freq >> ihh1 >> ihh2 >> data >> normedData >> critBool;

                while (pos > winEnd)
                {
                    winStarts[i].push_back(winStart);
                    nSNPs[i].push_back(numSNPsWin);
                    if (numSNPsWin == 0) fracCrit[i].push_back(-1);
                    else fracCrit[i].push_back(double(numCrit) / double(numSNPsWin));
                    if (numSNPsWin >= minSNPs && numCrit >= 0) numWindows++;
                    maxAbsScore[i].push_back(maxAbs);

                    maxAbs = -99999.9;
                    winStart += winSize;
                    winEnd += winSize;
                    numSNPsWin = 0;
                    numCrit = 0;
                }
                if (abs(normedData) > maxAbs) maxAbs = abs(normedData);
                numSNPsWin++;
                numCrit += critBool;
            }
            else if (statType == XPEHH)
            {
                fin >> name >> pos >> gpos >> freq1 >> ihh1 >> freq2 >> ihh2 >> data >> normedData >> critInt;

                while (pos > winEnd)
                {
                    winStarts[i].push_back(winStart);
                    nSNPs[i].push_back(numSNPsWin);

                    if (numSNPsWin < minSNPs)
                    {
                        fracCritTop[i].push_back(-1);
                        fracCritBot[i].push_back(-1);
                    }
                    else
                    {
                        fracCritTop[i].push_back(double(numCritTop) / double(numSNPsWin));
                        numWindowsTop++;
                        fracCritBot[i].push_back(double(numCritBot) / double(numSNPsWin));
                        numWindowsBot++;
                    }
                    maxScore[i].push_back(maxVal);
                    minScore[i].push_back(minVal);

                    maxVal = -99999.9;
                    minVal = 99999.9;
                    winStart += winSize;
                    winEnd += winSize;
                    numSNPsWin = 0;
                    numCritTop = 0;
                    numCritBot = 0;
                }
                if (normedData > maxVal) maxVal = normedData;
                if (normedData < minVal) minVal = normedData;
                numSNPsWin++;
                if (critInt == 1) numCritTop++;
                else if (critInt == -1) numCritBot++;
            }
            else if (statType == IHH12)
            {
                fin >> name >> pos >> freq1 >> data >> normedData >> critBool;

                while (pos > winEnd)
                {
                    winStarts[i].push_back(winStart);
                    nSNPs[i].push_back(numSNPsWin);
                    if (numSNPsWin == 0) fracCrit[i].push_back(-1);
                    else fracCrit[i].push_back(double(numCrit) / double(numSNPsWin));
                    if (numSNPsWin >= minSNPs && numCrit >= 0) numWindows++;
                    winStart += winSize;
                    winEnd += winSize;
                    numSNPsWin = 0;
                    numCrit = 0;
                }
                numSNPsWin++;
                numCrit += critBool;
            }
        }
        fin.close();
    }

    // For iHS and IHH12:
    if (statType == IHS || statType == IHH12)
    {
        cerr << numWindows << " nonzero windows.\n";
        flog << numWindows << " nonzero windows.\n";

        double *allSNPsPerWindow = new double[numWindows];
        double *allFracCritPerWindow = new double[numWindows];

        int k = 0;
        //Load all num SNPs per window into a single double vector to determine quantile boundaries across
        for (int i = 0; i < nfiles; i++)
        {
            for (int j = 0; j < nSNPs[i].size(); j++)
            {
                if (nSNPs[i][j] >= minSNPs && fracCrit[i][j] >= 0)
                {
                    allSNPsPerWindow[k] = nSNPs[i][j];
                    allFracCritPerWindow[k] = fracCrit[i][j];
                    k++;
                }
            }
        }

        //Sort allSNPsPerWindow and rearrange allFracCritPerWindow based on that sorting
        gsl_sort2(allSNPsPerWindow, 1, allFracCritPerWindow, 1, numWindows);

        double *quantileBound = new double[numQuantiles];

        //determine quantile boundaries
        for (int i = 0; i < numQuantiles; i++)
        {
            quantileBound[i] = gsl_stats_quantile_from_sorted_data(allSNPsPerWindow, 1, numWindows, double(i + 1) / double(numQuantiles));
        }


        /*
        *instead of splitting into a mini vector for each quantile bin, just pass a reference to the
        *start of the slice plus its size to gsl_stats_quantile_from_sorted_data
        *will need the number of snps per quantile bin
        */
        int b = 0; // quantile bin index
        int count = 0; //number of SNPs in quantile, not necessarily equal across quantiles because of ties
        int start = 0; // start index of the current quantile bin in allFracCritPerWindow
        map<string, double> *topWindowBoundary = new map<string, double>[numQuantiles];


        //cerr << "\nnSNPs 0.1 0.5 1.0 5.0\n";
        //flog << "\nnSNPs 0.1 0.5 1.0 5.0\n";

        cerr << "\nnSNPs 1.0 5.0\n";
        flog << "\nnSNPs 1.0 5.0\n";

        for (int i = 0; i < numWindows; i++)
        {
            if (allSNPsPerWindow[i] <= quantileBound[b])
            {
                count++;
            }
            else
            {
                gsl_sort(&(allFracCritPerWindow[start]), 1, count);

                topWindowBoundary[b]["1.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.990);
                topWindowBoundary[b]["5.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.950);

                cerr << quantileBound[b] << " "
                     << topWindowBoundary[b]["1.0"] << " "
                     << topWindowBoundary[b]["5.0"] << endl;

                flog << quantileBound[b] << " "
                     << topWindowBoundary[b]["1.0"] << " "
                     << topWindowBoundary[b]["5.0"] << endl;

                start = i;
                count = 0;
                b++;
            }
        }

        gsl_sort(&(allFracCritPerWindow[start]), 1, count);
        topWindowBoundary[b]["1.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.990);
        topWindowBoundary[b]["5.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.950);

        cerr << quantileBound[b] << " "
             << topWindowBoundary[b]["1.0"] << " "
             << topWindowBoundary[b]["5.0"] << "\n\n";

        flog << quantileBound[b] << " "
             << topWindowBoundary[b]["1.0"] << " "
             << topWindowBoundary[b]["5.0"] << "\n\n";

        delete[] allSNPsPerWindow;
        delete[] allFracCritPerWindow;

        for (int i = 0; i < nfiles; i++)
        {
            fout.open(winfilename[i].c_str());
            if (fout.fail())
            {
                cerr << "ERROR: " << winfilename[i] << " " << strerror(errno);
                throw -1;
            }
            cerr << "Creating window file " << winfilename[i] << endl;
            flog << "Creating window file " << winfilename[i] << endl;

            for (int j = 0; j < nSNPs[i].size(); j++)
            {
                if (nSNPs[i][j] < minSNPs || fracCrit[i][j] < 0)
                {
                    fout << winStarts[i][j] << "\t" << winStarts[i][j] + winSize << "\t"
                         << nSNPs[i][j] << "\t" << fracCrit[i][j] << "\t-1";
                    if (statType == IHS)
                    {
                        if (maxAbsScore[i][j] == -99999.9)
                            fout << "\tNA\n";
                        else
                            fout << "\t" << maxAbsScore[i][j] << endl;
                    }
                    else // IHH12 does not have maxAbsScore
                    {
                        fout << "\n";
                    }
                    continue;
                }

                double percentile = 100.0;
                for (b = 0; b < numQuantiles; b++)
                {
                    if (nSNPs[i][j] <= quantileBound[b])
                        break;
                }

                if (fracCrit[i][j] >= topWindowBoundary[b]["5.0"] && fracCrit[i][j] < topWindowBoundary[b]["1.0"])
                    percentile = 5.0;
                else if (fracCrit[i][j] >= topWindowBoundary[b]["1.0"])
                    percentile = 1.0;

                fout << winStarts[i][j] << "\t" << winStarts[i][j] + winSize << "\t"
                     << nSNPs[i][j] << "\t" << fracCrit[i][j] << "\t" << percentile;

                if (statType == IHS)
                {
                    if (maxAbsScore[i][j] == -99999.9)
                        fout << "\tNA\n";
                    else
                        fout << "\t" << maxAbsScore[i][j] << endl;
                }
                else
                {
                    fout << "\n";
                }
            }
            fout.close();
        }

        delete[] quantileBound;
        delete[] topWindowBoundary;
    }
    else if (statType == XPEHH)
    {
        cerr << numWindowsTop << " windows with nSNPs >= " << minSNPs << " (Top).\n";
        cerr << numWindowsBot << " windows with nSNPs >= " << minSNPs << " (Bot).\n";
        flog << numWindowsTop << " windows with nSNPs >= " << minSNPs << " (Top).\n";
        flog << numWindowsBot << " windows with nSNPs >= " << minSNPs << " (Bot).\n";

        double *allSNPsPerWindowTop = new double[numWindowsTop];
        double *allFracCritPerWindowTop = new double[numWindowsTop];
        double *allSNPsPerWindowBot = new double[numWindowsBot];
        double *allFracCritPerWindowBot = new double[numWindowsBot];

        int kTop = 0;
        int kBot = 0;

        for (int i = 0; i < nfiles; i++)
        {
            for (int j = 0; j < nSNPs[i].size(); j++)
            {
                if (nSNPs[i][j] >= minSNPs && fracCritTop[i][j] >= 0)
                {
                    allSNPsPerWindowTop[kTop] = nSNPs[i][j];
                    allFracCritPerWindowTop[kTop] = fracCritTop[i][j];
                    kTop++;
                }
                if (nSNPs[i][j] >= minSNPs && fracCritBot[i][j] >= 0)
                {
                    allSNPsPerWindowBot[kBot] = nSNPs[i][j];
                    allFracCritPerWindowBot[kBot] = fracCritBot[i][j];
                    kBot++;
                }
            }
        }

        gsl_sort2(allSNPsPerWindowTop, 1, allFracCritPerWindowTop, 1, numWindowsTop);
        gsl_sort2(allSNPsPerWindowBot, 1, allFracCritPerWindowBot, 1, numWindowsBot);

        double *quantileBoundTop = new double[numQuantiles];
        double *quantileBoundBot = new double[numQuantiles];

        for (int i = 0; i < numQuantiles; i++)
        {
            quantileBoundTop[i] = gsl_stats_quantile_from_sorted_data(allSNPsPerWindowTop, 1, numWindowsTop, double(i + 1) / double(numQuantiles));
            quantileBoundBot[i] = gsl_stats_quantile_from_sorted_data(allSNPsPerWindowBot, 1, numWindowsBot, double(i + 1) / double(numQuantiles));
        }

        // Process Top quantiles
        int bTop = 0, countTop = 0, startTop = 0;
        map<string, double> *topWindowBoundaryTop = new map<string, double>[numQuantiles];

        cerr << "\nHigh Scores\nnSNPs 1.0 5.0\n";
        flog << "\nHigh Scores\nnSNPs 1.0 5.0\n";

        for (int i = 0; i < numWindowsTop; i++)
        {
            if (allSNPsPerWindowTop[i] <= quantileBoundTop[bTop])
            {
                countTop++;
            }
            else
            {
                gsl_sort(&(allFracCritPerWindowTop[startTop]), 1, countTop);

                topWindowBoundaryTop[bTop]["1.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowTop[startTop]), 1, countTop, 0.990);
                topWindowBoundaryTop[bTop]["5.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowTop[startTop]), 1, countTop, 0.950);

                cerr << quantileBoundTop[bTop] << " "
                     << topWindowBoundaryTop[bTop]["1.0"] << " "
                     << topWindowBoundaryTop[bTop]["5.0"] << endl;

                flog << quantileBoundTop[bTop] << " "
                     << topWindowBoundaryTop[bTop]["1.0"] << " "
                     << topWindowBoundaryTop[bTop]["5.0"] << endl;

                startTop = i;
                countTop = 0;
                bTop++;
            }
        }
        gsl_sort(&(allFracCritPerWindowTop[startTop]), 1, countTop);
        topWindowBoundaryTop[bTop]["1.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowTop[startTop]), 1, countTop, 0.990);
        topWindowBoundaryTop[bTop]["5.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowTop[startTop]), 1, countTop, 0.950);

        cerr << quantileBoundTop[bTop] << " "
             << topWindowBoundaryTop[bTop]["1.0"] << " "
             << topWindowBoundaryTop[bTop]["5.0"] << "\n\n";

        flog << quantileBoundTop[bTop] << " "
             << topWindowBoundaryTop[bTop]["1.0"] << " "
             << topWindowBoundaryTop[bTop]["5.0"] << "\n\n";

        // Process Bot quantiles
        int bBot = 0, countBot = 0, startBot = 0;
        map<string, double> *topWindowBoundaryBot = new map<string, double>[numQuantiles];

        cerr << "\nLow Scores\nnSNPs 1.0 5.0\n";
        flog << "\nLow Scores\nnSNPs 1.0 5.0\n";

        for (int i = 0; i < numWindowsBot; i++)
        {
            if (allSNPsPerWindowBot[i] <= quantileBoundBot[bBot])
            {
                countBot++;
            }
            else
            {
                gsl_sort(&(allFracCritPerWindowBot[startBot]), 1, countBot);

                topWindowBoundaryBot[bBot]["1.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowBot[startBot]), 1, countBot, 0.990);
                topWindowBoundaryBot[bBot]["5.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowBot[startBot]), 1, countBot, 0.950);

                cerr << quantileBoundBot[bBot] << " "
                     << topWindowBoundaryBot[bBot]["1.0"] << " "
                     << topWindowBoundaryBot[bBot]["5.0"] << endl;

                flog << quantileBoundBot[bBot] << " "
                     << topWindowBoundaryBot[bBot]["1.0"] << " "
                     << topWindowBoundaryBot[bBot]["5.0"] << endl;

                startBot = i;
                countBot = 0;
                bBot++;
            }
        }
        gsl_sort(&(allFracCritPerWindowBot[startBot]), 1, countBot);
        topWindowBoundaryBot[bBot]["1.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowBot[startBot]), 1, countBot, 0.990);
        topWindowBoundaryBot[bBot]["5.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowBot[startBot]), 1, countBot, 0.950);

        cerr << quantileBoundBot[bBot] << " "
             << topWindowBoundaryBot[bBot]["1.0"] << " "
             << topWindowBoundaryBot[bBot]["5.0"] << "\n\n";

        flog << quantileBoundBot[bBot] << " "
             << topWindowBoundaryBot[bBot]["1.0"] << " "
             << topWindowBoundaryBot[bBot]["5.0"] << "\n\n";

        delete[] allSNPsPerWindowTop;
        delete[] allFracCritPerWindowTop;
        delete[] allSNPsPerWindowBot;
        delete[] allFracCritPerWindowBot;

        // Now output files for each input file
        for (int i = 0; i < nfiles; i++)
        {
            fout.open(winfilename[i].c_str());
            if (fout.fail())
            {
                cerr << "ERROR: " << winfilename[i] << " " << strerror(errno);
                throw -1;
            }
            cerr << "Creating window file " << winfilename[i] << endl;
            flog << "Creating window file " << winfilename[i] << endl;

            for (int j = 0; j < nSNPs[i].size(); j++)
            {
                if (nSNPs[i][j] < minSNPs || fracCritTop[i][j] < 0 || fracCritBot[i][j] < 0)
                {
                    fout << winStarts[i][j] << "\t" << winStarts[i][j] + winSize << "\t"
                         << nSNPs[i][j] << "\t" << fracCritTop[i][j] << "\t" << fracCritBot[i][j] << "\t-1\t-1\n";
                    continue;
                }

                double percentileTop = 100.0, percentileBot = 100.0;
                int bTopOut, bBotOut;
                for (bTopOut = 0; bTopOut < numQuantiles; bTopOut++)
                {
                    if (nSNPs[i][j] <= quantileBoundTop[bTopOut]) break;
                }
                for (bBotOut = 0; bBotOut < numQuantiles; bBotOut++)
                {
                    if (nSNPs[i][j] <= quantileBoundBot[bBotOut]) break;
                }

                if (fracCritTop[i][j] >= topWindowBoundaryTop[bTopOut]["5.0"] && fracCritTop[i][j] < topWindowBoundaryTop[bTopOut]["1.0"])
                    percentileTop = 5.0;
                else if (fracCritTop[i][j] >= topWindowBoundaryTop[bTopOut]["1.0"])
                    percentileTop = 1.0;

                if (fracCritBot[i][j] >= topWindowBoundaryBot[bBotOut]["5.0"] && fracCritBot[i][j] < topWindowBoundaryBot[bBotOut]["1.0"])
                    percentileBot = 5.0;
                else if (fracCritBot[i][j] >= topWindowBoundaryBot[bBotOut]["1.0"])
                    percentileBot = 1.0;

                fout << winStarts[i][j] << "\t" << winStarts[i][j] + winSize << "\t"
                     << nSNPs[i][j] << "\t" << fracCritTop[i][j] << "\t" << fracCritBot[i][j]
                     << "\t" << percentileTop << "\t" << percentileBot << "\t"
                     << maxScore[i][j] << "\t" << minScore[i][j] << endl;
            }
            fout.close();
        }

        delete[] quantileBoundTop;
        delete[] quantileBoundBot;
        delete[] topWindowBoundaryTop;
        delete[] topWindowBoundaryBot;
    }

    // Clean up
    if (statType == IHS || statType == IHH12)
    {
        delete[] fracCrit;
        delete[] maxAbsScore;
    }
    else if (statType == XPEHH)
    {
        delete[] fracCritTop;
        delete[] fracCritBot;
        delete[] maxScore;
        delete[] minScore;
    }
    delete[] winStarts;
    delete[] nSNPs;
    delete[] winfilename;
}


/*
void analyzeSNPWindows(string normedfiles[],int fileLoci[], int nfiles, int snpWinSize)
{
  cerr << "\nAnalyzing SNP windows:\n\n";
  vector<int>* winStarts = new vector<int>[nfiles];
  vector<int>* winEnds = new vector<int>[nfiles];
  vector<string>* startSNP = new vector<int>[nfiles];
  vector<string>* endSNP = new vector<int>[nfiles];
  vector<double>* fracCrit = new vector<double>[nfiles];
  ifstream fin;
  ofstream fout;
  string* winfilename = new string[nfiles];

  char str[10];
  sprintf(str,"%d",snpWinSize);

  string name;
  int pos;
  double freq, ihh1, ihh2, data, normedData;
  bool crit;
  int numWindows = 0;

  for (int i = 0; i < nfiles; i++)
    {
      fin.open(normedfiles[i].c_str());
      if(fin.fail())
    {
      cerr << "ERROR: " << normedfiles[i] << " " << strerror(errno);
      throw -1;
    }

      //generate winfile names
      winfilename[i] = normedfiles[i];
      winfilename[i] += ".";
      winfilename[i] += str;
      winfilename[i] += "snp.windows";

      //Load information into vectors for analysis
      int winStart = 1;
      int winEnd = winStart + winSize - 1;
      int snpsInWin = 0;
      int numCrit = 0;
      for(int j = 0; j < fileLoci[i]; j++)
    {
      fin >> name;
      fin >> pos;
      fin >> freq;
      fin >> ihh1;
      fin >> ihh2;
      fin >> data;
      fin >> normedData;
      fin >> crit;

      snpsInWin++;
      numCrit+=crit;
    }
    }
}
*/

// return number of columns in file
int countCols(ifstream &fin)
{
    int current_cols = 0;
    int currentPos = fin.tellg();

    string line;
    // read the rest of this line
    getline(fin, line);
    // read the next full line
    getline(fin, line);

    current_cols = countFields(line);

    fin.clear();
    // restore the previous position
    fin.seekg(currentPos);

    currentPos = fin.tellg();

    return current_cols;
}

// return number of columns to skip
int colsToSkip(ifstream &fin, int numCols)
{
    // determine the number of cols to skip 
    // (the number more than the number we care about: 6)
    int presentNumCols = countCols(fin);
    int numberColsToSkip = 0;
    string junk;

    if ( presentNumCols > numCols)
    {
        numberColsToSkip = presentNumCols - numCols;
    }

    return numberColsToSkip;
}


// skip the first numCols columns in the file
void skipCols(ifstream &fin, int numCols)
{
    string junk;
    
    for(int i=0; i<numCols; i++)   
    {
        fin >> junk;
    }
}


// Calculates mean and variance for each bin, then normalizes data in each bin
// Also normalizes the full data array so that we can calculate quantiles later
// freq[] is the derived allele frequency array
// data[] is the raw iHS or XPEHH score array
// nloci is the number of loci in the arrays
// mean[] is an array of size numBins to hold the mean for each bin
// variance[] is an array of size numBins to hold the variance for each bin
// n[] is an array of size numBins to hold the number of loci in each bin
// numBins is the number of frequency bins
// threshold[] is an array of size numBins holding the upper boundary for each bin
void getMeanVarBins(double freq[], double data[], int nloci, double mean[], double variance[], int n[], int numBins, double threshold[], string normLogInput)
{

    //initialize
    for (int b = 0; b < numBins; b++)
    {
        n[b] = 0;
        mean[b] = 0;
        variance[b] = 0;
    }


    if(normLogInput!=""){

    }else{
        //Calculate sum(x_i) stored in mean[b], and sum(x_i^2) stored in variance[b] in each frequency bin b
        for (int i = 0; i < nloci; i++)
        {
            if (data[i] == MISSING) continue;
            for (int b = 0; b < numBins; b++)
            {
                if (freq[i] < threshold[b])
                {
                    n[b]++;
                    mean[b] += data[i];
                    variance[b] += data[i] * data[i];
                    break;
                }
            }
        }
    }
    
    //Transform the sum(x_i) and sum(x_i^2) into mean and variances for each bin
    double temp;
    for (int b = 0; b < numBins; b++)
    {
        temp = ( variance[b] - (mean[b] * mean[b]) / (n[b]) ) / (n[b] - 1);
        variance[b] = temp;
        temp = mean[b] / n[b];
        mean[b] = temp;
    }

    //normalize the full data array
    //so that we can calculate quntiles later
    for (int i = 0; i < nloci; i++)
    {
        if (data[i] == MISSING) continue;
        for (int b = 0; b < numBins; b++)
        {
            if (freq[i] < threshold[b])
            {
                data[i] = (data[i] - mean[b]) / sqrt(variance[b]);
                break;
            }
        }
    }
    return;
}


/**
 * Unified normalization by allele-frequency bins for XP-EHH / XPNSL, IHH12, and IHS.
 *
 * type:
 *   - "xp"     : expects columns  name pos gpos freq1 ihh1 freq2 ihh2 data
 *   - "ihh12"  : expects columns  name pos freq1 data
 *   - "ihs"    : expects (at least) name pos freq ihh1 ihh2 ihs (+ optional extra cols)
 *
 * Behavior matches your original functions:
 *   - "xp"  writes header + "\tnormxpehh|normxpnsl\tcrit\n" (depending on XPNSL).
 *   - "ihh12" writes header + "\tnormihh12\tcrit\n".
 *   - "ihs"  determines extra columns via colsToSkip(fin, 6) and does NOT write a header.
 */
void normalizeDataByBins(const string& filename,
                     const string& outfilename,
                     const int& fileLoci,
                     double mean[], double variance[], int n[],
                     int numBins, double threshold[],
                     double upperCutoff, double lowerCutoff,
                     const string& type,
                     bool XPNSL)
{
    std::ifstream fin(filename.c_str());
    std::ofstream fout(outfilename.c_str());
    if (!fin) {
        throw std::runtime_error("ERROR opening input: " + filename + " : " + std::strerror(errno));
    }
    if (!fout) {
        throw std::runtime_error("ERROR opening output: " + outfilename + " : " + std::strerror(errno));
    }

    if (type != "xp" && type != "ihh12" && type != "ihs") {
        throw std::invalid_argument("normalizeByBins: type must be one of {\"xp\",\"ihh12\",\"ihs\"}");
    }

    string name, header, junk;
    int pos;
    double gpos, freq1, freq2, ihh1, ihh2, data, normedData;

    if (type == "xp") {
        // XP-EHH / XPNSL variant
        std::getline(fin, header);
        if (XPNSL) fout << header << "\t" << "normxpnsl" << "\t" << "crit" << "\n";
        else       fout << header << "\t" << "normxpehh" << "\t" << "crit" << "\n";

        for (int j = 0; j < fileLoci; ++j) {
            int numInBin = 0;  // reset per row
            fin >> name >> pos >> gpos >> freq1 >> ihh1 >> freq2 >> ihh2 >> data;
            if (!fin) break;

            if (data == MISSING) continue;

            // find bin
            for (int b = 0; b < numBins; ++b) {
                if (freq1 < threshold[b]) {
                    double v = variance[b];
                    if (v < 0) v = 0; // clamp tiny negatives
                    normedData = (v > 0) ? (data - mean[b]) / std::sqrt(v) : 0.0;
                    numInBin   = n[b];
                    break;
                }
            }

            if (numInBin >= 20) {
                fout << name << "\t"
                     << pos  << "\t"
                     << gpos << "\t"
                     << freq1 << "\t"
                     << ihh1 << "\t"
                     << freq2 << "\t"
                     << ihh2 << "\t"
                     << data  << "\t"
                     << normedData << "\t";
                if (normedData >= upperCutoff)      fout << "1\n";
                else if (normedData <= lowerCutoff) fout << "-1\n";
                else                                 fout << "0\n";
            }
        }

    } else if (type == "ihh12") {
        // IHH12 variant
        std::getline(fin, header);
        fout << header << "\t" << "normihh12" << "\t" << "crit" << "\n";

        for (int j = 0; j < fileLoci; ++j) {
            int numInBin = 0;  // reset per row
            fin >> name >> pos >> freq1 >> data;
            if (!fin) break;

            if (data == MISSING) continue;

            for (int b = 0; b < numBins; ++b) {
                if (freq1 < threshold[b]) {
                    double v = variance[b];
                    if (v < 0) v = 0;
                    normedData = (v > 0) ? (data - mean[b]) / std::sqrt(v) : 0.0;
                    numInBin   = n[b];
                    break;
                }
            }

            if (numInBin >= 20) {
                fout << name << "\t"
                     << pos  << "\t"
                     << freq1 << "\t"
                     << data  << "\t"
                     << normedData << "\t";
                if (normedData >= upperCutoff || normedData <= lowerCutoff) fout << "1\n";
                else                                                        fout << "0\n";
            }
        }

    } else { // type == "ihs"
        // IHS variant — determine extra columns beyond the 6 we care about
        int numColsToSkip = colsToSkip(fin, 6);

        for (int j = 0; j < fileLoci; ++j) {
            int numInBin = 0;  // reset per row
            double freq;
            fin >> name >> pos >> freq >> ihh1 >> ihh2 >> data;
            if (!fin) break;

            // skip any extra columns present
            if (numColsToSkip > 0) skipCols(fin, numColsToSkip);

            if (data == MISSING) continue;

            for (int b = 0; b < numBins; ++b) {
                if (freq < threshold[b]) {
                    double v = variance[b];
                    if (v < 0) v = 0;
                    normedData = (v > 0) ? (data - mean[b]) / std::sqrt(v) : 0.0;
                    numInBin   = n[b];
                    break;
                }
            }

            if (numInBin >= 20) {
                // Note: original IHS normalizer did not write a header; keep same output shape.
                fout << name << "\t"
                     << pos  << "\t"
                     << freq << "\t"
                     << ihh1 << "\t"
                     << ihh2 << "\t"
                     << data << "\t"
                     << normedData << "\t";
                if (normedData >= upperCutoff || normedData <= lowerCutoff) fout << "1\n";
                else                                                        fout << "0\n";
            }
        }
    }

    fin.close();
    fout.close();
}


// returns number of lines in file
// checks that each line has the expected number of columns
// if not, throws an error
int checkFile(std::ifstream &fin, const std::string &name)
{
    // Map file "name" to its allowed column counts
    // IHS supports an alternate format (e.g., with --ihs-detail)
    const static std::unordered_map<std::string, std::vector<int>> allowedCols = {
        {"ihh12", {4}},
        {"ihs",   {6, 10}},
        {"xpehh", {8}}
    };

    auto it = allowedCols.find(name);
    if (it == allowedCols.end())
    {
        std::cerr << "ERROR: unknown file type name \"" << name
                  << "\". Expected one of: IHH12, IHS, XPEHH.\n";
        throw 0;
    }
    const std::vector<int>& expected = it->second;

    std::string line;
    int current_cols = 0;

    // Save and restore stream position properly using streampos
    std::streampos start = fin.tellg();

    int nloci = 0;
    while (std::getline(fin, line))
    {
        nloci++;
        current_cols = countFields(line);

        // Skip header (first line) column check, as in your original code
        if (nloci > 1)
        {
            bool ok = std::find(expected.begin(), expected.end(), current_cols) != expected.end();
            if (!ok)
            {
                // Build a helpful "expected" message
                std::cerr << "ERROR: line " << nloci << " has " << current_cols
                          << " columns, but expected ";
                if (expected.size() == 1)
                {
                    std::cerr << expected[0] << " columns.\n";
                }
                else
                {
                    std::cerr << expected[0] << " or " << expected[1] << " columns.\n";
                }
                throw 0;
            }
        }
    }

    // Exclude header from count if file had at least one line
    if (nloci > 0) nloci--;

    // Reset the stream state and position
    fin.clear();
    fin.seekg(start);

    return nloci;
}




/**
 * Unified reader for IHH12, XP-EHH, and IHS tables.
 *
 * @param filenames  list of file paths
 * @param fileLoci   number of loci (rows) to read from each file
 * @param nfiles     number of files
 * @param type       one of: "ihh12", "xp", "ihs"
 * @param freq1      output array for freq (or pop1 freq)
 * @param score      output array for score
 * @param freq2      (optional) output array for pop2 freq; required for type="xp"
 */
void readAll(const vector<string>& filenames,
             const int fileLoci[],
             int nfiles,
             const string& type,
             double* freq1,
             double* score,
             double* freq2)
{
    if (type != "ihh12" && type != "xp" && type != "ihs") {
        throw std::invalid_argument("readAll: type must be one of {\"ihh12\",\"xp\",\"ihs\"}");
    }
    if (type == "xp" && freq2 == nullptr) {
        throw std::invalid_argument("readAll: freq2 must be non-null for type=\"xp\"");
    }

    ifstream fin;
    string junk;
    int overallCount = 0;

    for (int i = 0; i < nfiles; ++i) {
        fin.open(filenames[i].c_str());
        if (!fin) {
            throw std::runtime_error("readAll: failed to open file: " + filenames[i]);
        }

        if (type == "ihh12") {
            // Skip header line
            std::getline(fin, junk);

            // For each row: SNP, pos, freq1, score
            for (int j = 0; j < fileLoci[i]; ++j) {
                fin >> junk;                      // SNP
                fin >> junk;                      // pos
                fin >> freq1[overallCount];       // freq
                fin >> score[overallCount];       // ihh12 score
                ++overallCount;
            }
        } else if (type == "xp") {
            // Skip header line
            std::getline(fin, junk);

            // Expected order (matching your original reader):
            // SNP, pos, junk, freq1, junk, freq2, junk, score
            for (int j = 0; j < fileLoci[i]; ++j) {
                fin >> junk;                      // SNP
                fin >> junk;                      // pos
                fin >> junk;                      // (extra col)
                fin >> freq1[overallCount];       // pop1 freq
                fin >> junk;                      // (sep/label)
                fin >> freq2[overallCount];       // pop2 freq
                fin >> junk;                      // (sep/label)
                fin >> score[overallCount];       // XP-EHH
                ++overallCount;
            }
        } else { // type == "ihs"
            // Determine how many extra columns are present beyond the expected 6
            // (SNP, pos, freq, ihh1, ihh2, ihs). Extra cols will be skipped per row.
            int numColsToSkip = colsToSkip(fin, 6);

            for (int j = 0; j < fileLoci[i]; ++j) {
                fin >> junk;                      // SNP
                fin >> junk;                      // pos
                fin >> freq1[overallCount];       // freq
                fin >> junk;                      // ihh1
                fin >> junk;                      // ihh2
                fin >> score[overallCount];       // ihs
                if (numColsToSkip > 0) {
                    skipCols(fin, numColsToSkip);
                }
                ++overallCount;
            }
        }

        fin.close();
    }
}


int countFields(const string &str)
{
    string::const_iterator it;
    int result;
    int numFields = 0;
    int seenChar = 0;
    for (it = str.begin() ; it < str.end(); it++)
    {
        result = isspace(*it);
        if (result == 0 && seenChar == 0)
        {
            numFields++;
            seenChar = 1;
        }
        else if (result != 0)
        {
            seenChar = 0;
        }
    }
    return numFields;
}

bool isint(string str)
{
    for (string::iterator it = str.begin(); it != str.end(); it++)
    {
        if (!isdigit(*it)) return 0;
    }

    return 1;
}
