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
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include "param_t.h"

using namespace std;

const string VERSION = "1.3.0";

const string PREAMBLE = " -- a program for downstream analysis of selscan output\n\
Source code and binaries can be found at\n\
\t<https://www.github.com/szpiech/selscan>\n\
\n\
Citations:\n\
\n\
selscan: ZA Szpiech and RD Hernandez (2014) MBE, 31: 2824-2827.\n\
iHH12: R Torres, et al. (2017) bioRxiv, doi: https://doi.org/10.1101/181859.\n\
       N Garud, et al. (2015) PLoS Genetics, 11: 1–32.\n\
nSL: A Ferrer-Admetlla, et al. (2014) MBE, 31: 1275-1291.\n\
xpehh: PC Sabeti, et al. (2007) Nature, 449: 913–918.\n\
iHS: BF Voight, et al. (2006) PLoS Biology, 4: e72.\n\
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


const int MISSING = -9999;

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
void normalizeXPEHHDataByBins(string &filename, string &outfilename, int &fileLoci, double mean[], double variance[], int n[], int numBins, double threshold[], double upperCutoff, double lowerCutoff);
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

int main(int argc, char *argv[])
{
    cerr << "norm v" + VERSION + "\n";
    param_t params;
    params.setPreamble(PREAMBLE);
    params.addFlag(ARG_FREQ_BINS, DEFAULT_FREQ_BINS, "", HELP_FREQ_BINS);
    params.addListFlag(ARG_FILES, DEFAULT_FILES, "", HELP_FILES);
    params.addFlag(ARG_LOG, DEFAULT_LOG, "", HELP_LOG);
    params.addFlag(ARG_WINSIZE, DEFAULT_WINSIZE, "", HELP_WINSIZE);
    params.addFlag(ARG_QBINS, DEFAULT_QBINS, "", HELP_QBINS);
    params.addFlag(ARG_MINSNPS, DEFAULT_MINSNPS, "", HELP_MINSNPS);
    params.addFlag(ARG_SNPWIN, DEFAULT_SNPWIN, "SILENT", HELP_SNPWIN);
    params.addFlag(ARG_SNPWINSIZE, DEFAULT_SNPWINSIZE, "SILENT", HELP_SNPWINSIZE);
    params.addFlag(ARG_BPWIN, DEFAULT_BPWIN, "", HELP_BPWIN);
    params.addFlag(ARG_FIRST, DEFAULT_FIRST, "", HELP_FIRST);
    params.addFlag(ARG_CRIT_NUM, DEFAULT_CRIT_NUM, "", HELP_CRIT_NUM);
    params.addFlag(ARG_CRIT_PERCENT, DEFAULT_CRIT_PERCENT, "", HELP_CRIT_PERCENT);
    params.addFlag(ARG_IHS, DEFAULT_IHS, "", HELP_IHS);
    params.addFlag(ARG_NSL, DEFAULT_NSL, "", HELP_NSL);
    params.addFlag(ARG_SOFT, DEFAULT_SOFT, "", HELP_SOFT);
    params.addFlag(ARG_XPEHH, DEFAULT_XPEHH, "", HELP_XPEHH);
    params.addFlag(ARG_XPNSL, DEFAULT_XPNSL, "", HELP_XPNSL);


    try
    {
        params.parseCommandLine(argc, argv);
    }
    catch (...)
    {
        return 1;
    }

    int numBins = params.getIntFlag(ARG_FREQ_BINS);
    vector<string> filename = params.getStringListFlag(ARG_FILES);
    int nfiles = filename.size();
    int winSize = params.getIntFlag(ARG_WINSIZE);
    string infoOutfile = params.getStringFlag(ARG_LOG);
    int numQBins = params.getIntFlag(ARG_QBINS);
    int minSNPs = params.getIntFlag(ARG_MINSNPS);
    int snpWinSize = params.getIntFlag(ARG_SNPWINSIZE);
    bool BPWIN = params.getBoolFlag(ARG_BPWIN);
    bool SNPWIN = params.getBoolFlag(ARG_SNPWIN);
    bool FIRST = params.getBoolFlag(ARG_FIRST);
    double critNum = params.getDoubleFlag(ARG_CRIT_NUM);
    double critPercent = params.getDoubleFlag(ARG_CRIT_PERCENT);
    bool IHS = params.getBoolFlag(ARG_IHS);
    bool NSL = params.getBoolFlag(ARG_NSL);
    bool SOFT = params.getBoolFlag(ARG_SOFT);
    bool XPEHH = params.getBoolFlag(ARG_XPEHH);
    bool XPNSL = params.getBoolFlag(ARG_XPNSL);

    if(XPNSL) XPEHH = true;

    if (numBins <= 0)
    {
        cerr << "ERROR: Must have a positive integer of frequency bins.\n";
        return 1;
    }

    if (numQBins <= 0)
    {
        cerr << "ERROR: Must have a positive integer of quantile bins.\n";
        return 1;
    }

    if (winSize <= 0)
    {
        cerr << "ERROR: Must have a positive integer window size.\n";
        return 1;
    }

    if (critNum <= 0)
    {
        cerr << "ERROR: Must give a positive critical value for |iHS| scores.\n";
        return 1;
    }

    if (critPercent != DEFAULT_CRIT_PERCENT && (critPercent <= 0 || critPercent >= 1))
    {
        cerr << "ERROR: Critical percentile must be in (0,1).\n";
        return 1;
    }

    
    if(IHS + XPEHH + NSL + SOFT!= 1){
        cerr << "Must specify exactly one of " + ARG_IHS + ", " + ARG_XPEHH + "," + ARG_NSL + "," + ARG_SOFT + "," + ARG_XPNSL + ".\n";
        return 1;
    }
    cerr << "You have provided " << nfiles << " output files for joint normalization.\n";

    string *outfilename = new string[nfiles];
    int *fileLoci = new int[nfiles];

    //ifstream* fin = new ifstream[nfiles];
    //ofstream* fout = new ofstream[nfiles];

    ifstream fin;

    int totalLoci = 0;

    //logging
    flog.open(infoOutfile.c_str());
    if (flog.fail())
    {
        cerr << "ERROR: " << infoOutfile << " " << strerror(errno) << endl;
        return 1;
    }

    for (int i = 0; i < argc; i++)
    {
        flog << argv[i] << " ";
    }
    flog << "\n\n";

    //flog << "Input files:\n";

    //For each file, open it, and check it for integrity
    //Also record total number of lines so we can allocate
    //enough space for the array of paired data that will
    //be used to calculate E[X] and E[X^2]
    for (int i = 0; i < nfiles; i++)
    {
        char str[10];
        sprintf(str, "%d", numBins);
        if(IHS || NSL) outfilename[i] = filename[i] + "." + str + "bins.norm";
        if(XPEHH || SOFT) outfilename[i] = filename[i] + ".norm";
        
        fin.open(filename[i].c_str());
        if (fin.fail())
        {
            cerr << "ERROR: " << infoOutfile << " " << strerror(errno);
            flog << "ERROR: " << infoOutfile << " " << strerror(errno);
            return 1;
        }
        else
        {
            cerr << "Opened " << filename[i] << endl;
            //flog << filename[i] << endl;
        }

        //check integrity of file and keep count of the number of lines
        try
        {
            if (IHS || NSL) fileLoci[i] = checkIHSfile(fin);
            if (SOFT) fileLoci[i] = checkIHH12file(fin);
            if (XPEHH) fileLoci[i] = checkXPEHHfile(fin);
            totalLoci += fileLoci[i];
        }
        catch (...)
        {
            return 1;
        }
        fin.close();
    }

    cerr << "\nTotal loci: " << totalLoci << endl;
    flog << "\nTotal loci: " << totalLoci << endl;

    if (IHS || NSL)
    {
        cerr << "Reading all data.\n";
        double *freq = new double[totalLoci];
        double *score = new double[totalLoci];
        //read in all data
        readAllIHS(filename, fileLoci, nfiles, freq, score);

        double *mean = new double[numBins];
        double *variance = new double[numBins];
        int *n = new int[numBins];

        double minFreq;
        double maxFreq;
        double step;

        //This would use the empirical range to draw bin boundaries
        //gsl_stats_minmax(&minFreq,&maxFreq,freq,1,totalLoci);

        //This uses the possible range to draw bin boundaries
        minFreq = 0.0;
        maxFreq = 1.0;

        step = (maxFreq - minFreq) / double(numBins);

        double *threshold = new double[numBins];

        for (int b = 0; b < numBins; b++)
        {
            threshold[b] = minFreq + (b + 1) * step;
        }

        cerr << "Calculating mean and variance per frequency bin:\n\n";
        getMeanVarBins(freq, score, totalLoci, mean, variance, n, numBins, threshold);

        gsl_sort(score, 1, totalLoci);

        double upperCutoff, lowerCutoff;

        if (critPercent != DEFAULT_CRIT_PERCENT && (critPercent > 0 && critPercent < 1))
        {
            upperCutoff = gsl_stats_quantile_from_sorted_data (score, 1, totalLoci, 1 - critPercent / 2.0 );
            lowerCutoff = gsl_stats_quantile_from_sorted_data (score, 1, totalLoci, critPercent / 2.0);

            cerr << "\nTop cutoff: " << upperCutoff << endl;
            cerr << "Bottom cutoff: " << lowerCutoff << "\n\n";
            flog << "\nTop cutoff: " << upperCutoff << endl;
            flog << "Bottom cutoff: " << lowerCutoff << "\n\n";
        }
        else
        {
            upperCutoff = critNum;
            lowerCutoff = -critNum;
        }
        delete [] freq;
        delete [] score;

        //Output bins info to file.
        cerr << "bin\tnum\tmean\tvariance\n";
        flog << "bin\tnum\tmean\tvariance\n";
        for (int i = 0; i < numBins; i++)
        {
            cerr << threshold[i] << "\t" << n[i] <<  "\t" << mean[i] << "\t" << variance[i] << endl;
            flog << threshold[i] << "\t" << n[i] <<  "\t" << mean[i] << "\t" << variance[i] << endl;
        }


        //Read each file and create normed files.
        if (FIRST) nfiles = 1;
        for (int i = 0; i < nfiles; i++)
        {
            cerr << "Normalizing " << filename[i] << "\n";
            normalizeIHSDataByBins(filename[i], outfilename[i], fileLoci[i], mean, variance, n, numBins, threshold, upperCutoff, lowerCutoff);
            //fin[i].close();
            //fout[i].close();
        }

        delete [] threshold;
        delete [] mean;
        delete [] variance;
        delete [] n;

        if (BPWIN) analyzeIHSBPWindows(outfilename, fileLoci, nfiles, winSize, numQBins, minSNPs);
        //if(SNPWIN) analyzeSNPWindows(outfilename,fileLoci,nfiles,snpWinSize);
    }
    else if (XPEHH || SOFT) {

        cerr << "Reading all data.\n";
        double *freq1 = new double[totalLoci];
        double *freq2;
        if(XPEHH) freq2 = new double[totalLoci];
        double *score = new double[totalLoci];
        //read in all data
        if(XPEHH) readAllXPEHH(filename, fileLoci, nfiles, freq1, freq2, score);
        if(SOFT) readAllIHH12(filename, fileLoci, nfiles, freq1, score);
        numBins = 1;
        double *mean = new double[numBins];
        double *variance = new double[numBins];
        int *n = new int[numBins];

        double minFreq;
        double maxFreq;
        double step;

        //This would use the empirical range to draw bin boundaries
        //gsl_stats_minmax(&minFreq,&maxFreq,freq,1,totalLoci);

        //This uses the possible range to draw bin boundaries
        minFreq = 0.0;
        maxFreq = 1.0;

        step = (maxFreq - minFreq) / double(numBins);

        double *threshold = new double[numBins];

        for (int b = 0; b < numBins; b++)
        {
            threshold[b] = minFreq + (b + 1) * step;
        }

        cerr << "Calculating mean and variance:\n\n";
        getMeanVarBins(freq1, score, totalLoci, mean, variance, n, numBins, threshold);

        gsl_sort(score, 1, totalLoci);

        double upperCutoff, lowerCutoff;

        if (critPercent != DEFAULT_CRIT_PERCENT && (critPercent > 0 && critPercent < 1))
        {
            upperCutoff = gsl_stats_quantile_from_sorted_data (score, 1, totalLoci, 1 - critPercent / 2.0 );
            lowerCutoff = gsl_stats_quantile_from_sorted_data (score, 1, totalLoci, critPercent / 2.0);

            cerr << "\nTop cutoff: " << upperCutoff << endl;
            cerr << "Bottom cutoff: " << lowerCutoff << "\n\n";
            flog << "\nTop cutoff: " << upperCutoff << endl;
            flog << "Bottom cutoff: " << lowerCutoff << "\n\n";
        }
        else
        {
            upperCutoff = critNum;
            lowerCutoff = -critNum;
        }
        delete [] freq1;
        if (XPEHH) delete [] freq2;
        delete [] score;

        //Output bins info to file.
        cerr << "num\tmean\tvariance\n";
        flog << "num\tmean\tvariance\n";
        for (int i = 0; i < numBins; i++)
        {
            cerr << n[i] <<  "\t" << mean[i] << "\t" << variance[i] << endl;
            flog << n[i] <<  "\t" << mean[i] << "\t" << variance[i] << endl;
        }


        //Read each file and create normed files.
        if (FIRST) nfiles = 1;
        for (int i = 0; i < nfiles; i++)
        {
            cerr << "Normalizing " << filename[i] << "\n";
            if(XPEHH) normalizeXPEHHDataByBins(filename[i], outfilename[i], fileLoci[i], mean, variance, n, numBins, threshold, upperCutoff, lowerCutoff);
            if(SOFT) normalizeIHH12DataByBins(filename[i], outfilename[i], fileLoci[i], mean, variance, n, numBins, threshold, upperCutoff, lowerCutoff);
            //fin[i].close();
            //fout[i].close();
        }

        delete [] threshold;
        delete [] mean;
        delete [] variance;
        delete [] n;

        if(BPWIN){
            if (XPEHH) analyzeXPEHHBPWindows(outfilename, fileLoci, nfiles, winSize, numQBins, minSNPs);
            if (SOFT) analyzeIHH12BPWindows(outfilename, fileLoci, nfiles, winSize, numQBins, minSNPs);
        }
    }
    flog.close();
    return 0;
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

void skipCols(ifstream &fin, int numCols)
{
    string junk;
    
    for(int i=0; i<numCols; i++)   
    {
        fin >> junk;
    }
}

void analyzeIHSBPWindows(string normedfiles[], int fileLoci[], int nfiles, int winSize, int numQuantiles, int minSNPs)
{
    cerr << "\nAnalyzing BP windows:\n\n";
    //int totalLoci = 0;
    //for (int i = 0; i < nfiles; i++) totalLoci+=fileLoci[i];
    vector<int> *winStarts = new vector<int>[nfiles];
    vector<int> *nSNPs = new vector<int>[nfiles];
    vector<double> *fracCrit = new vector<double>[nfiles];
    vector<double> *maxAbsScore = new vector<double>[nfiles];

    ifstream fin;
    ofstream fout;
    string *winfilename = new string[nfiles];

    char str[10];
    sprintf(str, "%d", winSize / 1000);

    string name;
    int pos;
    double freq, ihh1, ihh2, data, normedData;
    bool crit;
    int numWindows = 0;

    for (int i = 0; i < nfiles; i++)
    {
        fin.open(normedfiles[i].c_str());
        if (fin.fail())
        {
            cerr << "ERROR: " << normedfiles[i] << " " << strerror(errno);
            throw - 1;
        }

        //generate winfile names
        winfilename[i] = normedfiles[i];
        winfilename[i] += ".";
        winfilename[i] += str;
        winfilename[i] += "kb.windows";

        //Load information into vectors for analysis
        int winStart = 1;
        int winEnd = winStart + winSize - 1;
        int numSNPs = 0;
        int numCrit = 0;
        int maxAbs = -99999;
        for (int j = 0; j < fileLoci[i]; j++)
        {
            fin >> name;
            fin >> pos;
            fin >> freq;
            fin >> ihh1;
            fin >> ihh2;
            fin >> data;
            fin >> normedData;
            fin >> crit;

            while (pos > winEnd)
            {
                winStarts[i].push_back(winStart);
                nSNPs[i].push_back(numSNPs);
                if (numSNPs == 0) fracCrit[i].push_back(-1);
                else fracCrit[i].push_back(double(numCrit) / double(numSNPs));

                if (numSNPs >= minSNPs && numCrit >= 0) numWindows++;
                maxAbsScore[i].push_back(maxAbs);
                
                maxAbs = -99999;
                winStart += winSize;
                winEnd += winSize;
                numSNPs = 0;
                numCrit = 0;
            }

            if(abs(normedData) > maxAbs) maxAbs = abs(normedData);
            numSNPs++;
            numCrit += crit;
        }
        fin.close();
    }

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
        quantileBound[i] = gsl_stats_quantile_from_sorted_data (allSNPsPerWindow, 1, numWindows, double(i + 1) / double(numQuantiles));
    }

    /*
     *instead of splitting into a mini vector for each quantile bin, just pass a reference to the
     *start of the slice plus its size to gsl_stats_quantile_from_sorted_data
     *will need the number of snps per quantile bin
     */
    int b = 0;//quantileBoundary index
    int count = 0;//number in quantile, not necessarily equal across quantiles because of ties
    int start = 0;//starting index for the sort function
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

            //topWindowBoundary[b]["0.1"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.999);
            //topWindowBoundary[b]["0.5"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.995);
            topWindowBoundary[b]["1.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.990);
            topWindowBoundary[b]["5.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.950);

            cerr << quantileBound[b] << " "
                 //<< topWindowBoundary[b]["0.1"] << " "
                 //<< topWindowBoundary[b]["0.5"] << " "
                 << topWindowBoundary[b]["1.0"] << " "
                 << topWindowBoundary[b]["5.0"] << endl;

            flog << quantileBound[b] << " "
                 //<< topWindowBoundary[b]["0.1"] << " "
                 //<< topWindowBoundary[b]["0.5"] << " "
                 << topWindowBoundary[b]["1.0"] << " "
                 << topWindowBoundary[b]["5.0"] << endl;

            start = i;
            count = 0;
            b++;
        }
    }

    gsl_sort(&(allFracCritPerWindow[start]), 1, count);
    //topWindowBoundary[b]["0.1"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.999);
    //topWindowBoundary[b]["0.5"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.995);
    topWindowBoundary[b]["1.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.990);
    topWindowBoundary[b]["5.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.950);

    cerr << quantileBound[b] << " "
         //<< topWindowBoundary[b]["0.1"] << " "
         //<< topWindowBoundary[b]["0.5"] << " "
         << topWindowBoundary[b]["1.0"] << " "
         << topWindowBoundary[b]["5.0"] << "\n\n";

    flog << quantileBound[b] << " "
         //<< topWindowBoundary[b]["0.1"] << " "
         //<< topWindowBoundary[b]["0.5"] << " "
         << topWindowBoundary[b]["1.0"] << " "
         << topWindowBoundary[b]["5.0"] << "\n\n";

    delete [] allSNPsPerWindow;
    delete [] allFracCritPerWindow;

    for (int i = 0; i < nfiles; i++)
    {
        fout.open(winfilename[i].c_str());
        if (fout.fail())
        {
            cerr << "ERROR: " << winfilename[i] << " " << strerror(errno);
            throw - 1;
        }
        cerr << "Creating window file " << winfilename[i] << endl;
        flog << "Creating window file " << winfilename[i] << endl;
        for (int j = 0; j < nSNPs[i].size(); j++)
        {
            if (nSNPs[i][j] < minSNPs || fracCrit[i][j] < 0)
            {
                fout << winStarts[i][j] << "\t" << winStarts[i][j] + winSize << "\t" << nSNPs[i][j] << "\t" << fracCrit[i][j] << "\t-1\tNA" << endl;
                continue;
            }
            double percentile = 100.0;
            for (b = 0; b < numQuantiles; b++)
            {
                if (nSNPs[i][j] <= quantileBound[b]) break;
            }

            if (fracCrit[i][j] >= topWindowBoundary[b]["5.0"] && fracCrit[i][j] < topWindowBoundary[b]["1.0"])
            {
                percentile = 5.0;
            }
            else if (fracCrit[i][j] >= topWindowBoundary[b]["1.0"])// && fracCrit[i][j] < topWindowBoundary[b]["0.5"])
            {
                percentile = 1.0;
            }
            /*
            else if (fracCrit[i][j] >= topWindowBoundary[b]["0.5"] && fracCrit[i][j] < topWindowBoundary[b]["0.1"])
            {
                percentile = 0.5;
            }
            else if (fracCrit[i][j] >= topWindowBoundary[b]["0.1"])
            {
                percentile = 0.1;
            }
            */
            fout << winStarts[i][j] << "\t" << winStarts[i][j] + winSize << "\t" << nSNPs[i][j] << "\t" << fracCrit[i][j] << "\t" << percentile << "\t";
            if(maxAbsScore[i][j] == -99999){
                fout << "NA\n";
            }
            else{
                fout << maxAbsScore[i][j] << endl;
            }
        }
        fout.close();
    }

    delete [] quantileBound;
    delete [] topWindowBoundary;
    delete [] winStarts;
    delete [] nSNPs;
    delete [] fracCrit;
    delete [] winfilename;

    return;
}

void analyzeXPEHHBPWindows(string normedfiles[], int fileLoci[], int nfiles, int winSize, int numQuantiles, int minSNPs)
{
    cerr << "\nAnalyzing BP windows:\n\n";
    //int totalLoci = 0;
    //for (int i = 0; i < nfiles; i++) totalLoci+=fileLoci[i];
    vector<int> *winStarts = new vector<int>[nfiles];
    vector<int> *nSNPs = new vector<int>[nfiles];
    vector<double> *fracCritTop = new vector<double>[nfiles];
    vector<double> *fracCritBot = new vector<double>[nfiles];
    vector<double> *maxScore = new vector<double>[nfiles];
    vector<double> *minScore = new vector<double>[nfiles];

    ifstream fin;
    ofstream fout;
    string *winfilename = new string[nfiles];

    char str[10];
    sprintf(str, "%d", winSize / 1000);

    string name, header;
    int pos;
    double gpos, freq1, freq2, ihh1, ihh2, data, normedData;
    int crit;
    int numWindowsTop = 0;
    int numWindowsBot = 0;

    for (int i = 0; i < nfiles; i++)
    {
        fin.open(normedfiles[i].c_str());
        if (fin.fail())
        {
            cerr << "ERROR: " << normedfiles[i] << " " << strerror(errno);
            throw - 1;
        }

        getline(fin, header);

        //generate winfile names
        winfilename[i] = normedfiles[i];
        winfilename[i] += ".";
        winfilename[i] += str;
        winfilename[i] += "kb.windows";

        //Load information into vectors for analysis
        int winStart = 1;
        int winEnd = winStart + winSize - 1;
        int numSNPs = 0;
        int numCritTop = 0;
        int numCritBot = 0;
        double max = -99999;
        double min = 99999;
        for (int j = 0; j < fileLoci[i]; j++)
        {
            fin >> name;
            fin >> pos;
            fin >> gpos;
            fin >> freq1;
            fin >> ihh1;
            fin >> freq2;
            fin >> ihh2;
            fin >> data;
            fin >> normedData;
            fin >> crit;

            while (pos > winEnd)
            {
                winStarts[i].push_back(winStart);
                nSNPs[i].push_back(numSNPs);
                if (numSNPs < minSNPs){
                    fracCritTop[i].push_back(-1);
                    fracCritBot[i].push_back(-1);
                }
                else{
                    fracCritTop[i].push_back(double(numCritTop) / double(numSNPs));
                    numWindowsTop++;
                    fracCritBot[i].push_back(double(numCritBot) / double(numSNPs));
                    numWindowsBot++;
                }
                maxScore[i].push_back(max);
                minScore[i].push_back(min);
                
                max = -99999;
                min = 99999;
                winStart += winSize;
                winEnd += winSize;
                numSNPs = 0;
                numCritTop = 0;
                numCritBot = 0;
            }

            if(normedData > max) max = normedData;
            if(normedData < min) min = normedData;
            numSNPs++;
            if(crit == 1) numCritTop++;
            else if (crit == -1) numCritBot++;
        }
        fin.close();
    }

    cerr << numWindowsTop << " windows with nSNPs >= " << minSNPs << ".\n";
    flog << numWindowsTop << " windows with nSNPs >= " << minSNPs << ".\n";
    double *allSNPsPerWindowTop = new double[numWindowsTop];
    double *allFracCritPerWindowTop = new double[numWindowsTop];

    double *allSNPsPerWindowBot = new double[numWindowsBot];
    double *allFracCritPerWindowBot = new double[numWindowsBot];

    int kTop = 0;
    int kBot = 0;
    //Load all num SNPs per window into a single double vector to determine quantile boundaries across
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

    //Sort allSNPsPerWindow and rearrange allFracCritPerWindow based on that sorting
    gsl_sort2(allSNPsPerWindowTop, 1, allFracCritPerWindowTop, 1, numWindowsTop);
    gsl_sort2(allSNPsPerWindowBot, 1, allFracCritPerWindowBot, 1, numWindowsBot);

    double *quantileBoundTop = new double[numQuantiles];
    double *quantileBoundBot = new double[numQuantiles];

    //determine quantile boundaries
    for (int i = 0; i < numQuantiles; i++)
    {
        quantileBoundTop[i] = gsl_stats_quantile_from_sorted_data (allSNPsPerWindowTop, 1, numWindowsTop, double(i + 1) / double(numQuantiles));
        quantileBoundBot[i] = gsl_stats_quantile_from_sorted_data (allSNPsPerWindowBot, 1, numWindowsBot, double(i + 1) / double(numQuantiles));
    }



////TOP
    /*
     *instead of splitting into a mini vector for each quantile bin, just pass a reference to the
     *start of the slice plus its size to gsl_stats_quantile_from_sorted_data
     *will need the number of snps per quantile bin
     */
    int bTop = 0;//quantileBoundary index
    int countTop = 0;//number in quantile, not necessarily equal across quantiles because of ties
    int startTop = 0;//starting index for the sort function
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

    delete [] allSNPsPerWindowTop;
    delete [] allFracCritPerWindowTop;

///BOT 
/*
     *instead of splitting into a mini vector for each quantile bin, just pass a reference to the
     *start of the slice plus its size to gsl_stats_quantile_from_sorted_data
     *will need the number of snps per quantile bin
     */
    int bBot = 0;//quantileBoundary index
    int countBot = 0;//number in quantile, not necessarily equal across quantiles because of ties
    int startBot = 0;//starting index for the sort function
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

    delete [] allSNPsPerWindowBot;
    delete [] allFracCritPerWindowBot;


    for (int i = 0; i < nfiles; i++)
    {
        fout.open(winfilename[i].c_str());
        if (fout.fail())
        {
            cerr << "ERROR: " << winfilename[i] << " " << strerror(errno);
            throw - 1;
        }
        cerr << "Creating window file " << winfilename[i] << endl;
        flog << "Creating window file " << winfilename[i] << endl;
        for (int j = 0; j < nSNPs[i].size(); j++)
        {
            fout << winStarts[i][j] << "\t" << winStarts[i][j] + winSize << "\t" << nSNPs[i][j] << "\t" << fracCritTop[i][j] << "\t" << fracCritBot[i][j] << "\t";
            if (nSNPs[i][j] < minSNPs)
            {
                fout << "-1\t-1\tNA\tNA" << endl;
                continue;
            }

            double percentile = 100.0;
            for (bTop = 0; bTop < numQuantiles; bTop++)
            {
                if (nSNPs[i][j] <= quantileBoundTop[bTop]) break;
            }

            if (fracCritTop[i][j] >= topWindowBoundaryTop[bTop]["5.0"] && fracCritTop[i][j] < topWindowBoundaryTop[bTop]["1.0"])
            {
                percentile = 5.0;
            }
            else if (fracCritTop[i][j] >= topWindowBoundaryTop[bTop]["1.0"])// && fracCritTop[i][j] < topWindowBoundaryTop[b]["0.5"])
            {
                percentile = 1.0;
            }
            
            fout << percentile << "\t";

            percentile = 100.0;
            for (bBot = 0; bBot < numQuantiles; bBot++)
            {
                if (nSNPs[i][j] <= quantileBoundBot[bBot]) break;
            }

            if (fracCritBot[i][j] >= topWindowBoundaryBot[bBot]["5.0"] && fracCritBot[i][j] < topWindowBoundaryBot[bBot]["1.0"])
            {
                percentile = 5.0;
            }
            else if (fracCritBot[i][j] >= topWindowBoundaryBot[bBot]["1.0"])// && fracCritTop[i][j] < topWindowBoundaryTop[b]["0.5"])
            {
                percentile = 1.0;
            }
            
            fout << percentile << "\t";

            if(maxScore[i][j] == -99999){
                fout << "NA\t";
            }
            else{
                fout << maxScore[i][j] << "\t";
            }
            if(minScore[i][j] == 99999){
                fout << "NA\n";
            }
            else{
                fout << minScore[i][j] << endl;
            }
        }
        fout.close();
    }

    delete [] quantileBoundTop;
    delete [] topWindowBoundaryTop;
    delete [] quantileBoundBot;
    delete [] topWindowBoundaryBot;

    delete [] winStarts;
    delete [] nSNPs;
    delete [] fracCritTop;
    delete [] fracCritBot;

    delete [] winfilename;

    return;
}

void analyzeIHH12BPWindows(string normedfiles[], int fileLoci[], int nfiles, int winSize, int numQuantiles, int minSNPs)
{
    cerr << "\nAnalyzing BP windows:\n\n";
    //int totalLoci = 0;
    //for (int i = 0; i < nfiles; i++) totalLoci+=fileLoci[i];
    vector<int> *winStarts = new vector<int>[nfiles];
    vector<int> *nSNPs = new vector<int>[nfiles];
    vector<double> *fracCrit = new vector<double>[nfiles];

    ifstream fin;
    ofstream fout;
    string *winfilename = new string[nfiles];

    char str[10];
    sprintf(str, "%d", winSize / 1000);

    string name, header;
    int pos;
    double gpos, freq1, ihh12, data, normedData;
    bool crit;
    int numWindows = 0;

    for (int i = 0; i < nfiles; i++)
    {
        fin.open(normedfiles[i].c_str());
        if (fin.fail())
        {
            cerr << "ERROR: " << normedfiles[i] << " " << strerror(errno);
            throw - 1;
        }

        getline(fin, header);

        //generate winfile names
        winfilename[i] = normedfiles[i];
        winfilename[i] += ".";
        winfilename[i] += str;
        winfilename[i] += "kb.windows";

        //Load information into vectors for analysis
        int winStart = 1;
        int winEnd = winStart + winSize - 1;
        int numSNPs = 0;
        int numCrit = 0;
        for (int j = 0; j < fileLoci[i]; j++)
        {
            fin >> name;
            fin >> pos;
            fin >> freq1;
            fin >> data;
            fin >> normedData;
            fin >> crit;

            while (pos > winEnd)
            {
                winStarts[i].push_back(winStart);
                nSNPs[i].push_back(numSNPs);
                if (numSNPs == 0) fracCrit[i].push_back(-1);
                else fracCrit[i].push_back(double(numCrit) / double(numSNPs));

                if (numSNPs >= minSNPs && numCrit >= 0) numWindows++;

                winStart += winSize;
                winEnd += winSize;
                numSNPs = 0;
                numCrit = 0;
            }

            numSNPs++;
            numCrit += crit;
        }
        fin.close();
    }

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
        quantileBound[i] = gsl_stats_quantile_from_sorted_data (allSNPsPerWindow, 1, numWindows, double(i + 1) / double(numQuantiles));
    }

    /*
     *instead of splitting into a mini vector for each quantile bin, just pass a reference to the
     *start of the slice plus its size to gsl_stats_quantile_from_sorted_data
     *will need the number of snps per quantile bin
     */
    int b = 0;//quantileBoundary index
    int count = 0;//number in quantile, not necessarily equal across quantiles because of ties
    int start = 0;//starting index for the sort function
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

            //topWindowBoundary[b]["0.1"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.999);
            //topWindowBoundary[b]["0.5"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.995);
            topWindowBoundary[b]["1.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.990);
            topWindowBoundary[b]["5.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.950);

            cerr << quantileBound[b] << " "
                 //<< topWindowBoundary[b]["0.1"] << " "
                 //<< topWindowBoundary[b]["0.5"] << " "
                 << topWindowBoundary[b]["1.0"] << " "
                 << topWindowBoundary[b]["5.0"] << endl;

            flog << quantileBound[b] << " "
                 //<< topWindowBoundary[b]["0.1"] << " "
                 //<< topWindowBoundary[b]["0.5"] << " "
                 << topWindowBoundary[b]["1.0"] << " "
                 << topWindowBoundary[b]["5.0"] << endl;

            start = i;
            count = 0;
            b++;
        }
    }

    gsl_sort(&(allFracCritPerWindow[start]), 1, count);
    //topWindowBoundary[b]["0.1"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.999);
    //topWindowBoundary[b]["0.5"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.995);
    topWindowBoundary[b]["1.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.990);
    topWindowBoundary[b]["5.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.950);

    cerr << quantileBound[b] << " "
         //<< topWindowBoundary[b]["0.1"] << " "
         //<< topWindowBoundary[b]["0.5"] << " "
         << topWindowBoundary[b]["1.0"] << " "
         << topWindowBoundary[b]["5.0"] << "\n\n";

    flog << quantileBound[b] << " "
         //<< topWindowBoundary[b]["0.1"] << " "
         //<< topWindowBoundary[b]["0.5"] << " "
         << topWindowBoundary[b]["1.0"] << " "
         << topWindowBoundary[b]["5.0"] << "\n\n";

    delete [] allSNPsPerWindow;
    delete [] allFracCritPerWindow;

    for (int i = 0; i < nfiles; i++)
    {
        fout.open(winfilename[i].c_str());
        if (fout.fail())
        {
            cerr << "ERROR: " << winfilename[i] << " " << strerror(errno);
            throw - 1;
        }
        cerr << "Creating window file " << winfilename[i] << endl;
        flog << "Creating window file " << winfilename[i] << endl;
        for (int j = 0; j < nSNPs[i].size(); j++)
        {
            if (nSNPs[i][j] < minSNPs || fracCrit[i][j] < 0)
            {
                fout << winStarts[i][j] << "\t" << winStarts[i][j] + winSize << "\t" << nSNPs[i][j] << "\t" << fracCrit[i][j] << "\t-1" << endl;
                continue;
            }
            double percentile = 100.0;
            for (b = 0; b < numQuantiles; b++)
            {
                if (nSNPs[i][j] <= quantileBound[b]) break;
            }

            if (fracCrit[i][j] >= topWindowBoundary[b]["5.0"] && fracCrit[i][j] < topWindowBoundary[b]["1.0"])
            {
                percentile = 5.0;
            }
            else if (fracCrit[i][j] >= topWindowBoundary[b]["1.0"])// && fracCrit[i][j] < topWindowBoundary[b]["0.5"])
            {
                percentile = 1.0;
            }
            /*
            else if (fracCrit[i][j] >= topWindowBoundary[b]["0.5"] && fracCrit[i][j] < topWindowBoundary[b]["0.1"])
            {
                percentile = 0.5;
            }
            else if (fracCrit[i][j] >= topWindowBoundary[b]["0.1"])
            {
                percentile = 0.1;
            }
            */
            fout << winStarts[i][j] << "\t" << winStarts[i][j] + winSize << "\t" << nSNPs[i][j] << "\t" << fracCrit[i][j] << "\t" << percentile << endl;
        }
        fout.close();
    }

    delete [] quantileBound;
    delete [] topWindowBoundary;
    delete [] winStarts;
    delete [] nSNPs;
    delete [] fracCrit;
    delete [] winfilename;

    return;
}


void getMeanVarBins(double freq[], double data[], int nloci, double mean[], double variance[], int n[], int numBins, double threshold[])
{
    //initialize
    for (int b = 0; b < numBins; b++)
    {
        n[b] = 0;
        mean[b] = 0;
        variance[b] = 0;
    }

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

//Reads a file, calculates the normalized score, and
//outputs the original row plus normed score
void normalizeIHSDataByBins(string &filename, string &outfilename, int &fileLoci, double mean[], double variance[], int n[], int numBins, double threshold[], double upperCutoff, double lowerCutoff)
{
    ifstream fin;
    ofstream fout;

    fin.open(filename.c_str());
    fout.open(outfilename.c_str());
    if (fout.fail())
    {
        cerr << "ERROR: " << outfilename << " " << strerror(errno);
        throw 1;
    }

    string name;
    int pos;
    double freq, data, normedData, ihh1, ihh2;;
    int numInBin = 0;
    string junk;

    // determine the number of cols to skip 
    // (the number more than the number we care about: 6)
    int numColsToSkip = 0;
    numColsToSkip = colsToSkip(fin, 6);

    for (int j = 0; j < fileLoci; j++)
    {
        fin >> name;
        fin >> pos;
        fin >> freq;
        fin >> ihh1;
        fin >> ihh2;
        fin >> data;

        // read in and skip extra columns
        skipCols(fin, numColsToSkip);

        if (data == MISSING) continue;
        for (int b = 0; b < numBins; b++)
        {
            if (freq < threshold[b])
            {
                normedData = (data - mean[b]) / sqrt(variance[b]);
                numInBin = n[b];
                break;
            }
        }

        if (numInBin >= 20)
        {
            fout << name << "\t"
                 << pos << "\t"
                 << freq << "\t"
                 << ihh1 << "\t"
                 << ihh2 << "\t"
                 << data << "\t"
                 << normedData << "\t";
            if (normedData >= upperCutoff || normedData <= lowerCutoff) fout << "1\n";
            else fout << "0\n";
        }
    }

    fin.close();
    fout.close();

    return;
}

void normalizeXPEHHDataByBins(string &filename, string &outfilename, int &fileLoci, double mean[], double variance[], int n[], int numBins, double threshold[], double upperCutoff, double lowerCutoff)
{
    ifstream fin;
    ofstream fout;

    fin.open(filename.c_str());
    fout.open(outfilename.c_str());
    if (fout.fail())
    {
        cerr << "ERROR: " << outfilename << " " << strerror(errno);
        throw 1;
    }

    string name, header;
    int pos;
    double gpos, freq1, freq2, data, normedData, ihh1, ihh2;;
    int numInBin = 0;

    getline(fin, header);

    fout << header + "\tnormxpehh\tcrit\n";

    for (int j = 0; j < fileLoci; j++)
    {
        fin >> name;
        fin >> pos;
        fin >> gpos;
        fin >> freq1;
        fin >> ihh1;
        fin >> freq2;
        fin >> ihh2;
        fin >> data;

        if (data == MISSING) continue;
        for (int b = 0; b < numBins; b++)
        {
            if (freq1 < threshold[b])
            {
                normedData = (data - mean[b]) / sqrt(variance[b]);
                numInBin = n[b];
                break;
            }
        }

        if (numInBin >= 20)
        {
            fout << name << "\t"
                 << pos << "\t"
                 << gpos << "\t"
                 << freq1 << "\t"
                 << ihh1 << "\t"
                 << freq2 << "\t"
                 << ihh2 << "\t"
                 << data << "\t"
                 << normedData << "\t";
            if (normedData >= upperCutoff) fout << "1\n";
            else if (normedData <= lowerCutoff) fout << "-1\n";
            else fout << "0\n";
        }
    }

    fin.close();
    fout.close();

    return;
}

void normalizeIHH12DataByBins(string &filename, string &outfilename, int &fileLoci, double mean[], double variance[], int n[], int numBins, double threshold[], double upperCutoff, double lowerCutoff)
{
    ifstream fin;
    ofstream fout;

    fin.open(filename.c_str());
    fout.open(outfilename.c_str());
    if (fout.fail())
    {
        cerr << "ERROR: " << outfilename << " " << strerror(errno);
        throw 1;
    }

    string name, header;
    int pos;
    double gpos, freq1, data, normedData, ihh1, ihh2;;
    int numInBin = 0;

    getline(fin, header);

    fout << header + "\tnormihh12\tcrit\n";

    for (int j = 0; j < fileLoci; j++)
    {
        fin >> name;
        fin >> pos;
        fin >> freq1;
        fin >> data;

        if (data == MISSING) continue;
        for (int b = 0; b < numBins; b++)
        {
            if (freq1 < threshold[b])
            {
                normedData = (data - mean[b]) / sqrt(variance[b]);
                numInBin = n[b];
                break;
            }
        }

        if (numInBin >= 20)
        {
            fout << name << "\t"
                 << pos << "\t"
                 << freq1 << "\t"
                 << data << "\t"
                 << normedData << "\t";
            if (normedData >= upperCutoff || normedData <= lowerCutoff) fout << "1\n";
            else fout << "0\n";
        }
    }

    fin.close();
    fout.close();

    return;
}

//returns number of lines in file
int checkIHSfile(ifstream &fin)
{
    string line;
    int expected_cols = 6;
    int expected_cols_alternate = 10; // this is the case if --ihs-detail is specified (four extra columns for iHH left/right and ancestral/derived)
    int current_cols = 0;

    //beginning of the file stream
    int start = fin.tellg();

    int nloci = 0;
    while (getline(fin, line))
    {
        nloci++;
        current_cols = countFields(line);
        if ((current_cols != expected_cols && current_cols != expected_cols_alternate) && nloci > 1)
        {
            cerr << "ERROR: line " << nloci << " has " << current_cols
                 << " columns, but expected " << expected_cols << " or " << expected_cols_alternate << " columns.\n";
            throw 0;
        }
        //previous_cols = current_cols;
    }

    fin.clear();
    fin.seekg(start);

    return nloci;
}

void readAllIHS(vector<string> filename, int fileLoci[], int nfiles, double freq[], double score[])
{
    ifstream fin;
    string junk;
    int overallCount = 0;
    for (int i = 0; i < nfiles; i++)
    {
        fin.open(filename[i].c_str());

        int numColsToSkip = 0;
        numColsToSkip = colsToSkip(fin, 6);

        for (int j = 0; j < fileLoci[i]; j++)
        {
            fin >> junk;
            fin >> junk;
            fin >> freq[overallCount];
            fin >> junk;
            fin >> junk;
            fin >> score[overallCount];
            skipCols(fin, numColsToSkip);
            overallCount++;
        }
        fin.close();
    }

    return;
}

int checkXPEHHfile(ifstream &fin)
{
    string line;
    int expected_cols = 8;
    int current_cols = 0;

    //beginning of the file stream
    int start = fin.tellg();

    int nloci = 0;
    while (getline(fin, line))
    {
        nloci++;
        current_cols = countFields(line);
        if ((current_cols != expected_cols) && nloci > 1)
        {
            cerr << "ERROR: line " << nloci << " has " << current_cols
                 << " columns, but expected " << expected_cols << " columns.\n";
            throw 0;
        }
        //previous_cols = current_cols;
    }

    nloci--;

    fin.clear();
    fin.seekg(start);

    return nloci;
}


int checkIHH12file(ifstream &fin)
{
    string line;
    int expected_cols = 4;
    int current_cols = 0;

    //beginning of the file stream
    int start = fin.tellg();

    int nloci = 0;
    while (getline(fin, line))
    {
        nloci++;
        current_cols = countFields(line);
        if ((current_cols != expected_cols) && nloci > 1)
        {
            cerr << "ERROR: line " << nloci << " has " << current_cols
                 << " columns, but expected " << expected_cols << " columns.\n";
            throw 0;
        }
        //previous_cols = current_cols;
    }

    nloci--;

    fin.clear();
    fin.seekg(start);

    return nloci;
}

void readAllXPEHH(vector<string> filename, int fileLoci[], int nfiles, double freq1[], double freq2[], double score[])
{
    ifstream fin;
    string junk;
    int overallCount = 0;
    for (int i = 0; i < nfiles; i++)
    {
        fin.open(filename[i].c_str());
        getline(fin, junk);
        for (int j = 0; j < fileLoci[i]; j++)
        {
            fin >> junk;
            fin >> junk;
            fin >> junk;
            fin >> freq1[overallCount];
            fin >> junk;
            fin >> freq2[overallCount];
            fin >> junk;
            fin >> score[overallCount];
            overallCount++;
        }
        fin.close();
    }

    return;
}

void readAllIHH12(vector<string> filename, int fileLoci[], int nfiles, double freq1[], double score[])
{
    ifstream fin;
    string junk;
    int overallCount = 0;
    for (int i = 0; i < nfiles; i++)
    {
        fin.open(filename[i].c_str());
        getline(fin, junk);
        for (int j = 0; j < fileLoci[i]; j++)
        {
            fin >> junk;
            fin >> junk;
            fin >> freq1[overallCount];
            fin >> score[overallCount];
            overallCount++;
        }
        fin.close();
    }

    return;
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
