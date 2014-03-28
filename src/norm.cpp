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

const string PREAMBLE = "\nnorm -- a program for downstream analysis of selscan output\n\
Source code and binaries can be found at\n\
\t<https://www.github.com/szpiech/selscan>\n\
\n\
norm currently analyzes iHS output only at this time.\n\
\n\
Citations:\n\
\n\
ZA Szpiech and RD Hernandez (2014) arXiv:1403.6854 [q-bio.PE]\n\
BF Voight, et al. (2006) PLoS Biology, 4: e72.\n\
\n\
To normalize iHS output across frequency bins:\n\
\n\
./norm --files <file1.ihs.out> ... <fileN.ihs.out>\n\
\n\
To normalize iHS output and analyze non-overlapping windows of fixed bp for \n\
extreme iHS scores:\n\
\n\
./norm --files <file1.ihs.out> ... <fileN.ihs.out> --bp-win\n";

string ARG_FREQ_BINS = "--bins";
int DEFAULT_FREQ_BINS = 100;
string HELP_FREQ_BINS = "The number of frequency bins in [0,1] for score normalization.";

string ARG_FILES = "--files";
string DEFAULT_FILES = "infile";
string HELP_FILES = "A list of files delimited by whitespace for\n\
\tjoint normalization.\n\
\tExpected format: <locus name> <physical pos> <freq> <ihh1> <ihh2> <ihs>";

string ARG_LOG = "--log";
string DEFAULT_LOG = "logfile";
string HELP_LOG = "The log file name.";

string ARG_WINSIZE = "--winsize";
int DEFAULT_WINSIZE = 100000;
string HELP_WINSIZE = "The non-overlapping window size for calculating the percentage\n\
\tof extreme SNPs.";

string ARG_QBINS = "--qbins";
int DEFAULT_QBINS = 20;
string HELP_QBINS = "Outlying windows are binned by number of sites within each\n\
\twindow.  This is the number of quantile bins to use.";

string ARG_MINSNPS = "--min-snps";
int DEFAULT_MINSNPS = 10;
string HELP_MINSNPS = "Only consider a bp window if it has at least this many SNPs.";

string ARG_SNPWIN = "--snp-win";
bool DEFAULT_SNPWIN = false;
string HELP_SNPWIN = "<not implemented> If set, will use windows of a constant\n\
\tSNP size with varying bp length.";

string ARG_SNPWINSIZE = "--snp-win-size";
int DEFAULT_SNPWINSIZE = 50;
string HELP_SNPWINSIZE = "<not implemented> The number of SNPs in a window.";

string ARG_BPWIN = "--bp-win";
bool DEFAULT_BPWIN = false;
string HELP_BPWIN = "If set, will use windows of a constant bp size with varying\n\
\tnumber of SNPs.";

string ARG_IHS = "--ihs";
bool DEFAULT_IHS = false;
string HELP_IHS = "Do iHS normalization.";

string ARG_SOFT = "--soft";
bool DEFAULT_SOFT = false;
string HELP_SOFT = "Do soft-iHS normalization.";

string ARG_FIRST = "--first";
bool DEFAULT_FIRST = false;
string HELP_FIRST = "Output only the first file's normalization.";

const int MISSING = -9999;

//returns number of lines in file
//throws 0 if the file fails
int checkIHSfile(ifstream &fin);

void readAll(vector<string> filename, int fileLoci[], int nfiles, double freq[], double score[]);

void getMeanVarBins(double freq[], double data[], int nloci, double mean[], double variance[], int n[], int numBins, double threshold[]);
void normalizeDataByBins(string &filename, string &outfilename, int &fileLoci, double mean[], double variance[], int n[], int numBins, double threshold[], double upperCutoff, double lowerCutoff);

void analyzeBPWindows(string normedfiles[], int fileLoci[], int nfiles, int winSize, int numQuantiles, int minSNPs);

int countFields(const string &str);
bool isint(string str);

ofstream flog;

int main(int argc, char *argv[])
{
    param_t params;
    params.setPreamble(PREAMBLE);
    params.addFlag(ARG_FREQ_BINS, DEFAULT_FREQ_BINS, "", HELP_FREQ_BINS);
    params.addListFlag(ARG_FILES, DEFAULT_FILES, "", HELP_FILES);
    params.addFlag(ARG_LOG, DEFAULT_LOG, "", HELP_LOG);
    params.addFlag(ARG_WINSIZE, DEFAULT_WINSIZE, "", HELP_WINSIZE);
    params.addFlag(ARG_QBINS, DEFAULT_QBINS, "", HELP_QBINS);
    params.addFlag(ARG_MINSNPS, DEFAULT_MINSNPS, "", HELP_MINSNPS);
    params.addFlag(ARG_SNPWIN, DEFAULT_SNPWIN, "", HELP_SNPWIN);
    params.addFlag(ARG_SNPWINSIZE, DEFAULT_SNPWINSIZE, "", HELP_SNPWINSIZE);
    params.addFlag(ARG_BPWIN, DEFAULT_BPWIN, "", HELP_BPWIN);
    params.addFlag(ARG_FIRST, DEFAULT_FIRST, "", HELP_FIRST);


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

    flog << "Input files:\n";

    //For each file, open it, and check it for integrity
    //Also record total number of lines so we can allocate
    //enough space for the array of paired data that will
    //be used to calculate E[X] and E[X^2]
    for (int i = 0; i < nfiles; i++)
    {
        char str[10];
        sprintf(str, "%d", numBins);
        outfilename[i] = filename[i] + "." + str + "bins.norm";

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
            flog << filename[i] << endl;
        }

        //check integrity of file and keep count of the number of lines
        try
        {
            fileLoci[i] = checkIHSfile(fin);
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

    flog << "\nOutput files:\n";
    for (int i = 0; i < nfiles; i++)
    {
        flog << outfilename[i] << endl;
    }
    flog << endl;
    bool NORM_IHS = true;
    if (NORM_IHS)
    {
        cerr << "Reading all frequency and iHS data.\n";
        double *freq = new double[totalLoci];
        double *score = new double[totalLoci];
        //read in all data
        readAll(filename, fileLoci, nfiles, freq, score);

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
        /*
        gsl_sort(score, 1, totalLoci);
        double upperCutoff = gsl_stats_quantile_from_sorted_data (score, 1, totalLoci, 0.995);
        double lowerCutoff = gsl_stats_quantile_from_sorted_data (score, 1, totalLoci, 0.005);

        cerr << "\nTop 0.5% cutoff: " << upperCutoff << endl;
        cerr << "Bottom 0.5% cutoff: " << lowerCutoff << "\n\n";
        flog << "\nTop 0.5% cutoff: " << upperCutoff << endl;
        flog << "Bottom 0.5% cutoff: " << lowerCutoff << "\n\n";
        */
        double upperCutoff = 2;
        double lowerCutoff = -2;
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
            normalizeDataByBins(filename[i], outfilename[i], fileLoci[i], mean, variance, n, numBins, threshold, upperCutoff, lowerCutoff);
            //fin[i].close();
            //fout[i].close();
        }

        delete [] threshold;
        delete [] mean;
        delete [] variance;
        delete [] n;

        if (BPWIN) analyzeBPWindows(outfilename, fileLoci, nfiles, winSize, numQBins, minSNPs);
        //if(SNPWIN) analyzeSNPWindows(outfilename,fileLoci,nfiles,snpWinSize);
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
void analyzeBPWindows(string normedfiles[], int fileLoci[], int nfiles, int winSize, int numQuantiles, int minSNPs)
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

                if (numSNPs >= minSNPs && numCrit > 0) numWindows++;

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
            if (nSNPs[i][j] >= minSNPs && fracCrit[i][j] > 0)
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
    cerr << "\nnSNPs 0.1 0.5 1.0 5.0\n";
    flog << "\nnSNPs 0.1 0.5 1.0 5.0\n";
    for (int i = 0; i < numWindows; i++)
    {
        if (allSNPsPerWindow[i] <= quantileBound[b])
        {
            count++;
        }
        else
        {
            gsl_sort(&(allFracCritPerWindow[start]), 1, count);

            topWindowBoundary[b]["0.1"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.999);
            topWindowBoundary[b]["0.5"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.995);
            topWindowBoundary[b]["1.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.990);
            topWindowBoundary[b]["5.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.950);

            cerr << quantileBound[b] << " "
                 << topWindowBoundary[b]["0.1"] << " "
                 << topWindowBoundary[b]["0.5"] << " "
                 << topWindowBoundary[b]["1.0"] << " "
                 << topWindowBoundary[b]["5.0"] << endl;

            flog << quantileBound[b] << " "
                 << topWindowBoundary[b]["0.1"] << " "
                 << topWindowBoundary[b]["0.5"] << " "
                 << topWindowBoundary[b]["1.0"] << " "
                 << topWindowBoundary[b]["5.0"] << endl;

            start = i;
            count = 0;
            b++;
        }
    }

    gsl_sort(&(allFracCritPerWindow[start]), 1, count);
    topWindowBoundary[b]["0.1"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.999);
    topWindowBoundary[b]["0.5"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.995);
    topWindowBoundary[b]["1.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.990);
    topWindowBoundary[b]["5.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.950);

    cerr << quantileBound[b] << " "
         << topWindowBoundary[b]["0.1"] << " "
         << topWindowBoundary[b]["0.5"] << " "
         << topWindowBoundary[b]["1.0"] << " "
         << topWindowBoundary[b]["5.0"] << "\n\n";

    flog << quantileBound[b] << " "
         << topWindowBoundary[b]["0.1"] << " "
         << topWindowBoundary[b]["0.5"] << " "
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
            if (nSNPs[i][j] < minSNPs || fracCrit[i][j] <= 0)
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
            else if (fracCrit[i][j] >= topWindowBoundary[b]["1.0"] && fracCrit[i][j] < topWindowBoundary[b]["0.5"])
            {
                percentile = 1.0;
            }
            else if (fracCrit[i][j] >= topWindowBoundary[b]["0.5"] && fracCrit[i][j] < topWindowBoundary[b]["0.1"])
            {
                percentile = 0.5;
            }
            else if (fracCrit[i][j] >= topWindowBoundary[b]["0.1"])
            {
                percentile = 0.1;
            }

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
void normalizeDataByBins(string &filename, string &outfilename, int &fileLoci, double mean[], double variance[], int n[], int numBins, double threshold[], double upperCutoff, double lowerCutoff)
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

    for (int j = 0; j < fileLoci; j++)
    {
        fin >> name;
        fin >> pos;
        fin >> freq;
        fin >> ihh1;
        fin >> ihh2;
        fin >> data;

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


//returns number of lines in file
int checkIHSfile(ifstream &fin)
{
    string line;
    int expected_cols = 6;
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

    fin.clear();
    fin.seekg(start);

    return nloci;
}

void readAll(vector<string> filename, int fileLoci[], int nfiles, double freq[], double score[])
{
    ifstream fin;
    string junk;
    int overallCount = 0;
    for (int i = 0; i < nfiles; i++)
    {
        fin.open(filename[i].c_str());

        for (int j = 0; j < fileLoci[i]; j++)
        {
            fin >> junk;
            fin >> junk;
            fin >> freq[overallCount];
            fin >> junk;
            fin >> junk;
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
