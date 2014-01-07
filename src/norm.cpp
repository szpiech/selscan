/*
 * Normalize iHS scores across multiple files (i.e. chromosomes)
 * Simple CLI: ./norm <bins> <outfile1> ... <outfile22>
 * First:
 *  Read in all data needed only for calculating E[X] and E[X^2] of the whole data
 * Second:
 *  Load each file separately and normalize, output results by simply adding a column with the normed score
 */


#include <iostream>
#include <fstream>
#include <map>
#include <cmath>
#include <cstdlib>
#include <cctype>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>

using namespace std;

const int MISSING = -9999;

struct data_t
{
  int *pos;
  double *freq;
  double *data;
  double *ihh1;
  double *ihh2;
  int nloci;
  int ncols;
  string *lociNames;
};

//returns number of lines in file
//throws 0 if the file fails
int checkIHSfile(ifstream &fin);

//void normalizeScores(double data[], int nloci);
//void readData(ifstream &fin, data_t &data);
void readAll(ifstream fin[],int fileLoci[], int nfiles, double freq[], double score[]);

void getMeanVarBins(double freq[], double data[], int nloci, double mean[], double variance[], int n[], int numBins, double threshold[]);
void normalizeDataByBins(ifstream &fin, ofstream &fout, int &fileLoci, double mean[], double variance[], int n[], int numBins, double threshold[],double upperCutoff,double lowerCutoff);

//void normalizeBins(double freq[], double data[], int nloci, int numBins);
int countFields(const string &str);
bool isint(string str);

int main(int argc, char* argv[])
{
  if (argc <= 2)
    {
      cerr << "USAGE: " << argv[0] << " <bins> <infile1> ... <infileN>\n";
      cerr << "Input files are expected to have the following format:\n"
	   << "\t<locus name> <positon> <freq> <ihh1> <ihh2> <ihs>\n";
      cerr << "Normalized scores are output to <infile1>.<bins>bins.norm ... <infileN>.<bins>bins.norm in the following format:\n"
	   << "\t<locus name> <positon> <freq> <ihh1> <ihh2> <ihs> <normed ihs>\n";
      return 0;
    }

  int nfiles = argc-2;
  cerr << "You have provided "<< nfiles << " iHS output files for joint normalization.\n";

  string* filename = new string[nfiles];
  string* outfilename = new string[nfiles];
  int* fileLoci = new int[nfiles];
  ifstream* fin = new ifstream[nfiles];
  ofstream* fout = new ofstream[nfiles];
  ofstream fout2;
  int totalLoci = 0;

  if(!isint(argv[1]))
    {
      cerr << "ERROR: Number of frequency bins must be an integer.\n";
      return 1;
    }

  int numBins = atoi(argv[1]);
  
  if(numBins <= 0)
    {
      cerr << "ERROR: Must have a positive integer of frequency bins.\n";
      return 1;
    }
  
  //For each file, open it, and check it for integrity
  //Also record total number of lines so we can allocate
  //enough space for the array of paired data that will
  //be used to calculate E[X] and E[X^2]
  for(int i = 0; i < nfiles; i++)
    {
      filename[i] = argv[i+2];
      fin[i].open(filename[i].c_str());
      if(fin[i].fail())
	{
	  cerr << "ERROR: Could not open " << filename << " for reading.\n";
	  return 1;
	}
      else
	{
	  cerr << "Opened " << filename[i] << endl;
	}

      //check integrity of file and keep count of the number of lines
      try
	{
	  fileLoci[i] = checkIHSfile(fin[i]);
	  totalLoci += fileLoci[i];
	}
      catch (...)
	{
	  return 1;
	}

    }



  string infoOutfile = filename[0];
  infoOutfile += ".bins.info";
  
  fout2.open(infoOutfile.c_str());
  
  if(fout2.fail())
    {
      cerr << "ERROR: Could not open " << infoOutfile << " for reading.\n";
      return 1;
    }
  fout2 << "Input files:\n";
  for(int i = 0; i < nfiles; i++)
    {
      fout2 << filename[i] << endl;
    }

  cerr << "Total loci: " << totalLoci << endl;
  fout2 << "Total loci: " << totalLoci << endl;


  for(int i = 0; i < nfiles; i++)
    {
      char str[10];
      sprintf(str,"%d",numBins);
      outfilename[i] = filename[i] + "." + str + "bins.norm";

      fout[i].open(outfilename[i].c_str());

      if(fout[i].fail())
	{
	  cerr << "ERROR: Could not open " << outfilename[i] << " for reading.\n";
	  return 1;
	}
    }

  fout2 << "\nOutput files:\n";
  for(int i = 0; i < nfiles; i++)
    {
      fout2 << outfilename[i] << endl;
    }
  fout2 << endl;

  cerr << "Reading all frequency and iHS data.\n";
  double* freq = new double[totalLoci];
  double* score = new double[totalLoci];
  //read in all data 
  readAll(fin,fileLoci,nfiles,freq,score);

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

  step = (maxFreq - minFreq)/double(numBins);
  
  double *threshold = new double[numBins];
  
  for (int b = 0; b < numBins; b++)
    {
      threshold[b] = minFreq + (b+1)*step;
    }

  cerr << "Calculating mean and variance per frequency bin.\n";
  getMeanVarBins(freq,score,totalLoci,mean,variance,n,numBins,threshold);

  gsl_sort(score, 1, totalLoci);
  double upperCutoff = gsl_stats_quantile_from_sorted_data (score, 1, totalLoci, 0.995);
  double lowerCutoff = gsl_stats_quantile_from_sorted_data (score, 1, totalLoci, 0.005);

  cerr << "Top 0.5% cutoff: " << upperCutoff << endl;
  cerr << "Bottom 0.5% cutoff: " << lowerCutoff << endl;
  fout2 << "Top 0.5% cutoff: " << upperCutoff << endl;
  fout2 << "Bottom 0.5% cutoff: " << lowerCutoff << endl;

  delete [] freq;
  delete [] score;
  
  //Output bins info to file.
  fout2 << "bin\tnum\tmean\tvariance\n";
  for(int i = 0; i < numBins; i++)
    {
      fout2 << threshold[i] << "\t" << n[i] <<  "\t" << mean[i] << "\t" << variance[i] << endl;
    }
  fout2.close();


  //Read each file and create normed files.
  for (int i = 0; i < nfiles; i++)
    {
      cerr << "Normalizing " << filename[i] << "\n";
      normalizeDataByBins(fin[i],fout[i],fileLoci[i],mean,variance,n,numBins,threshold,upperCutoff,lowerCutoff);
    }

  cerr << "Finished.\n";

  delete [] threshold;

  return 0;
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
      if(data[i] == MISSING) continue;
      for (int b = 0; b < numBins; b++)
	{
	  if(freq[i] < threshold[b])
	    {
	      n[b]++;
	      mean[b] += data[i];
	      variance[b] += data[i]*data[i];
	      break;
	    }
	}
    }
  
  //Transform the sum(x_i) and sum(x_i^2) into mean and variances for each bin
  double temp;
  for (int b = 0; b < numBins; b++)
    {
      temp = ( variance[b] - (mean[b]*mean[b])/(n[b]) ) / (n[b] - 1);
      variance[b] = temp;
      temp = mean[b]/n[b];
      mean[b] = temp;
    }

  //normalize the full data array
  //so that we can calculate quntiles later
  for (int i = 0; i < nloci; i++)
    {
      if(data[i] == MISSING) continue;
      for (int b = 0; b < numBins; b++)
	{
	  if(freq[i] < threshold[b])
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
void normalizeDataByBins(ifstream &fin, ofstream &fout, int &fileLoci, double mean[], double variance[], int n[], int numBins, double threshold[],double upperCutoff, double lowerCutoff)
{

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
      
      if(data == MISSING) continue;
      for (int b = 0; b < numBins; b++)
	{
	  if(freq < threshold[b])
	    {
	      normedData = (data - mean[b]) / sqrt(variance[b]);
	      numInBin = n[b];
	      break;
	    }
	}
      
      if(numInBin >= 20)
	{
	  fout << name << "\t" 
	       << pos << "\t" 
	       << freq << "\t" 
	       << ihh1 << "\t" 
	       << ihh2 << "\t" 
	       << data << "\t" 
	       << normedData << "\t";
	  if(normedData >= upperCutoff || normedData <= lowerCutoff) fout << "1\n";
	  else fout << "0\n";
	}
    }

  return;
}

/*
void normalizeBins(double freq[], double data[], int nloci, int numBins)
{
  
  double *mean = new double[numBins];
  double *variance = new double[numBins];
  int *n = new int[numBins];
  
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
      if(data[i] == MISSING) continue;
      for (int b = 0; b < numBins; b++)
	{
	  if(freq[i] < double(b+1)/double(numBins))
	    {
	      n[b]++;
	      mean[b] += data[i];
	      variance[b] += data[i]*data[i];
	      break;
	    }
	}
    }

  //Transform the sum(x_i) and sum(x_i^2) into mean and variances for each bin
  double temp;
  for (int b = 0; b < numBins; b++)
    {
      temp = ( variance[b] - (mean[b]*mean[b])/(n[b]) ) / (n[b] - 1);
      variance[b] = temp;
      temp = mean[b]/n[b];
      mean[b] = temp;
    }


  //For each locus, determine the frequency bin it belongs to
  //and transform the value
  for (int i = 0; i < nloci; i++)
    {
      if(data[i] == MISSING) continue;
      for (int b = 0; b < numBins; b++)
	{
	  if(freq[i] < double(b+1)/double(numBins))
	    {
	      temp = (data[i] - mean[b]) / sqrt(variance[b]);
	      data[i] = temp;
	      break;
	    }
	}
    }

  delete [] mean;
  delete [] variance;
  delete [] n;

}
*/

//returns number of lines in file
int checkIHSfile(ifstream &fin)
{
  string line;
  int expected_cols = 6;
  int current_cols = 0;

  //beginning of the file stream
  int start = fin.tellg();

  int nloci = 0;
  while(getline(fin,line))
    {
      nloci++;
      current_cols = countFields(line);
      if((current_cols != expected_cols) && nloci > 1)
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

void readAll(ifstream fin[], int fileLoci[], int nfiles, double freq[], double score[])
{
  string junk;
  int overallCount = 0;
  for (int i = 0; i < nfiles; i++)
    {
      int start = fin[i].tellg();

      for (int j = 0; j < fileLoci[i]; j++)
	{
	  fin[i] >> junk;
	  fin[i] >> junk;
	  fin[i] >> freq[overallCount];
	  fin[i] >> junk;
	  fin[i] >> junk;
	  fin[i] >> score[overallCount];
	  overallCount++;
	}
      
      fin[i].clear();
      fin[i].seekg(start);
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
      if(result == 0 && seenChar == 0)
	{
	  numFields++;
	  seenChar = 1;
	}
      else if(result != 0)
	{
	  seenChar = 0;
	}
    }
  return numFields;
}

bool isint(string str)
{
  for (string::iterator it=str.begin(); it!=str.end(); it++)
    {
      if(!isdigit(*it)) return 0;
    }

  return 1;
}
