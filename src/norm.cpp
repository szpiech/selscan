#include "norm.h"
#include <deque>

// #include <gsl/gsl_statistics_double.h>
// #include <gsl/gsl_interp.h>


//int main(int argc, char *argv[])
int SelscanNorm::runToolNorm(int argc, char *argv[])
{
    cerr << "selscan v" + VERSION + "\n";
    cerr << "Subcommand: norm\n";
    cerr << "================\n";
    param_t params;
    params.setPreamble(PREAMBLE);
    params.addFlag(ARG_FREQ_BINS, DEFAULT_FREQ_BINS, "", HELP_FREQ_BINS);
    params.addListFlag(ARG_FILES, DEFAULT_FILES, "", HELP_FILES);
    params.addFlag(ARG_LOG, DEFAULT_LOG, "", HELP_LOG);
    params.addFlag(ARG_WINSIZE, DEFAULT_WINSIZE, "", HELP_WINSIZE);
    params.addFlag(ARG_QBINS, DEFAULT_QBINS, "", HELP_QBINS);
    params.addFlag(ARG_MINSNPS, DEFAULT_MINSNPS, "", HELP_MINSNPS);
    // params.addFlag(ARG_SNPWIN, DEFAULT_SNPWIN, "SILENT", HELP_SNPWIN);
    // params.addFlag(ARG_SNPWINSIZE, DEFAULT_SNPWINSIZE, "SILENT", HELP_SNPWINSIZE);
    params.addFlag(ARG_BPWIN, DEFAULT_BPWIN, "", HELP_BPWIN);
    params.addFlag(ARG_FIRST, DEFAULT_FIRST, "", HELP_FIRST);
    params.addFlag(ARG_CRIT_NUM, DEFAULT_CRIT_NUM, "", HELP_CRIT_NUM);
    params.addFlag(ARG_CRIT_PERCENT, DEFAULT_CRIT_PERCENT, "", HELP_CRIT_PERCENT);
    params.addFlag(ARG_IHS, DEFAULT_IHS, "", HELP_IHS);
    params.addFlag(ARG_NSL, DEFAULT_NSL, "", HELP_NSL);
    params.addFlag(ARG_SOFT, DEFAULT_SOFT, "", HELP_SOFT);
    params.addFlag(ARG_XPEHH, DEFAULT_XPEHH, "", HELP_XPEHH);
    params.addFlag(ARG_XPNSL, DEFAULT_XPNSL, "", HELP_XPNSL);


    params.addFlag(ARG_LOG_INPUT, DEFAULT_LOG_INPUT, "", HELP_LOG_INPUT); //added in v3
    params.addFlag(ARG_FINE_PERCENTILE, DEFAULT_FINE_PERCENTILE, "", HELP_FINE_PERCENTILE); //added in v3
    params.addFlag(ARG_BED, DEFAULT_BED, "", HELP_BED); //added in v3
    params.addFlag(ARG_WIN_FILE, DEFAULT_WIN_FILE, "", HELP_WIN_FILE); //added in v3
    params.addFlag(ARG_GENE_SETA, DEFAULT_GENE_SETA, "", HELP_GENE_SETA); //added in v3
    params.addFlag(ARG_GENE_SETB, DEFAULT_GENE_SETB, "", HELP_GENE_SETB); //added in v3
    // params.addFlag(ARG_NO_HEADER, DEFAULT_NO_HEADER, "", HELP_NO_HEADER); //added in v3

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
    // int snpWinSize = params.getIntFlag(ARG_SNPWINSIZE);
        // bool SNPWIN = params.getBoolFlag(ARG_SNPWIN);
    bool BPWIN = params.getBoolFlag(ARG_BPWIN);

    bool FIRST = params.getBoolFlag(ARG_FIRST);
    double critNum = params.getDoubleFlag(ARG_CRIT_NUM);
    double critPercent = params.getDoubleFlag(ARG_CRIT_PERCENT);
    
    bool IHS = params.getBoolFlag(ARG_IHS);
    bool NSL = params.getBoolFlag(ARG_NSL);
    bool SOFT = params.getBoolFlag(ARG_SOFT);
    bool XPEHH = params.getBoolFlag(ARG_XPEHH);
    bool XPNSL = params.getBoolFlag(ARG_XPNSL);

    this->FINE_PERCENTILE = params.getBoolFlag(ARG_FINE_PERCENTILE);
    this->GENE_BED = params.getStringFlag(ARG_BED);
    this->USE_GENE_BED = params.flagWasSet(ARG_BED);


    // read window file
    vector<string> windowFile = params.getStringListFlag(ARG_WIN_FILE);
    //string windowFile = params.getStringFlag(ARG_WIN_FILE);
    string geneBedFile = params.getStringFlag(ARG_BED);
    bool ANNOTATE_WINDOWS = params.flagWasSet(ARG_WIN_FILE);// && params.flagWasSet(ARG_BED);
    string geneSetA = params.getStringFlag(ARG_GENE_SETA);
    string geneSetB = params.getStringFlag(ARG_GENE_SETB);
    bool PERMUTE_TEST = params.flagWasSet(ARG_GENE_SETA) && params.flagWasSet(ARG_GENE_SETB);
    

    ///// ***** BLOCK START : ALLOWING GENE BASED ANALYSIS ***** /////
    //
    //
    bool DO_NORM = (nfiles > 0);                       // if input stat files provided
    bool DO_WINDOW = (winSize > 0);                    // user gave --winsize → do window-based stats
    bool DO_GENE = this->USE_GENE_BED;                // user gave --bed
    bool HAVE_WINFILE = (ANNOTATE_WINDOWS);           // user gave --win-file
    //
    // CASE 1: Just normalization
    //
    if (DO_NORM && !DO_WINDOW && !DO_GENE) {
        cerr << "[selscan-norm] Running normalization only...\n";
        // return runNormalization(params);
    }

    //
    // CASE 2: Normalization + window
    //
    if (DO_NORM && DO_WINDOW && !DO_GENE) {
        cerr << "[selscan-norm] Running normalization + window...\n";
        // string winOut = runNormalization(params);
        // runWindow(winOut, winSize);
        return 0;
    }

    //
    // CASE 3: Already-normalized file, window only
    //
    if (!DO_NORM && DO_WINDOW && !DO_GENE) {
        cerr << "[selscan-norm] Running windowing only...\n";
        // runWindow(params.getStringListFlag(ARG_FILES)[0], winSize);
        return 0;
    }

    //
    // CASE 4: Normalization + window + gene
    //
    if (DO_NORM && DO_WINDOW && DO_GENE && !HAVE_WINFILE) {
        cerr << "[selscan-norm] Running normalization + window + gene...\n";
        // string winOut = runNormalization(params);
        // string winFile = runWindow(winOut, winSize);
        //annotateWindows(geneBedFile, {winFile}, (XPEHH || XPNSL));
        return 0;
    }

    //
    // CASE 5: Already-normalized + window + gene
    //
    if (!DO_NORM && DO_WINDOW && DO_GENE && !HAVE_WINFILE) {
        cerr << "[selscan-norm] Running window + gene...\n";
        // string winFile = runWindow(params.getStringListFlag(ARG_FILES)[0], winSize);
        // annotateWindows(geneBedFile, {winFile}, (XPEHH || XPNSL));
        return 0;
    }

    //
    // CASE 6: Gene annotation with precomputed window file
    //
    if (DO_GENE && HAVE_WINFILE && !DO_WINDOW) {
        cerr << "[selscan-norm] Gene annotation on precomputed windows...\n";
        // annotateWindows(geneBedFile, windowFile, (XPEHH || XPNSL));
        return 0;
    }



    if(ANNOTATE_WINDOWS){
        if(IHS + XPEHH + NSL + SOFT + XPNSL != 1){
            cerr << "ERROR: Must specify exactly one of " + ARG_IHS + ", " + ARG_XPEHH + "," + ARG_NSL + "," + ARG_SOFT + "," + ARG_XPNSL + ".\n";
            return 1;
        }
        
        if(!this->USE_GENE_BED){
            cerr << "ERROR: --gene-bed must be provided with --annotate-win to annotate windows.\n";
            return 1;
        }
        if(geneBedFile == DEFAULT_BED || ANNOTATE_WINDOWS == false || windowFile[0] == DEFAULT_WIN_FILE){
            cerr << "ERROR: --gene-bed and --annotate-win must be provided to annotate windows.\n";
            return 1;
        }
        cerr << "Annotating " << windowFile.size() <<" windows " << " with genes from " << geneBedFile << "\n";
        // string statName;
        // if(IHS) statName = "ihs";
        // if(NSL) statName = "nsl";
        // if(XPEHH) statName = "xpehh";
        // if(XPNSL) statName = "xpnsl";
        // if(SOFT) statName = "ihh12";

        vector<string> filename = params.getStringListFlag(ARG_FILES);
        int nfiles = filename.size();

        bool XP = (XPEHH || XPNSL);
        annotateWindows(geneBedFile, windowFile, XP);
        return 0;
    }

    if(PERMUTE_TEST){
        if(geneSetA == DEFAULT_GENE_SETA || geneSetB == DEFAULT_GENE_SETB){
            cerr << "ERROR: --gene-target and --gene-background must be provided to run permutation test.\n";
            return 1;
        }
        cerr << "Running permutation test with target genes from " << geneSetA << " and background genes from " << geneSetB << "\n";
        perm_test(geneSetA, geneSetB);
        return 0;
    }

    if(params.flagWasSet(ARG_FREQ_BINS) && params.flagWasSet(ARG_LOG_INPUT)){
        cerr << "ERROR: Options --log-input and --bins cannot be used together. "
         << "--log-input already provides frequency bin information.\n";
        exit(1);
    }

    if(params.flagWasSet(ARG_FREQ_BINS) && (params.flagWasSet(ARG_XPEHH) || params.flagWasSet(ARG_XPNSL) || params.flagWasSet(ARG_SOFT))){
        cerr << "ERROR: Options --bins and --xpehh/--xpnsl/--ihh12 cannot be used together. "
         << "Frequency bins are only used for iHS and nSL normalization.\n";
        exit(1);
    }
       
    if(params.flagWasSet(ARG_QBINS) && params.flagWasSet(ARG_LOG_INPUT)){
        cerr << "ERROR: Options --log-input and --qbins cannot be used together. "
         << "--log-input already provides quantile bin information.\n";
        exit(1);
    }

    bool LOG_INPUT =  params.flagWasSet(ARG_LOG_INPUT);



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
    
    if(IHS + XPEHH + NSL + SOFT + XPNSL != 1){
        cerr << "ERROR: Must specify exactly one of " + ARG_IHS + ", " + ARG_XPEHH + "," + ARG_NSL + "," + ARG_SOFT + "," + ARG_XPNSL + ".\n";
        return 1;
    }

    // cout<<("=== Parameters ===\n");
    // cout<<params.flagWasSet(ARG_FILES)<<endl;
    
    if( filename[0] == DEFAULT_FILES ){
        if( !ANNOTATE_WINDOWS ){
            cerr << "ERROR: Must provide at least one input file with " + ARG_FILES + " for normalization.\n";
            return 1;
        }
    }

    if(ANNOTATE_WINDOWS && params.flagWasSet(ARG_FILES)){
        cerr << "ERROR: Cannot provide input files with " + ARG_FILES + " when annotating windows.\n";
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

    cerr << "\nTotal loci: " << totalLoci << endl; // across all files
    flog << "\nTotal loci: " << totalLoci << endl; // across all files


    string chr;
    //@NORMLOGINPUT
    if (IHS || NSL)
    {
        cerr << "Reading all data.\n";
        double *freq = new double[totalLoci];
        double *score = new double[totalLoci];
        //read in all data
        readAllIHS(filename, fileLoci, nfiles, freq, score, chr);

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
        if(LOG_INPUT){
            cerr << "Using frequency bins from " << params.getStringFlag(ARG_LOG_INPUT) << "\n";
            flog << "Using frequency bins from " << params.getStringFlag(ARG_LOG_INPUT) << "\n";

            getMeanVarBinsFromLog(params.getStringFlag(ARG_LOG_INPUT),
                                       freq, score, totalLoci,
                                       mean, variance,  n, numBins, threshold, (XPNSL|| XPEHH|| SOFT));
        }else{
            getMeanVarBins(freq, score, totalLoci, mean, variance, n, numBins, threshold);

        }

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
            normalizeIHSDataByBins(filename[i], outfilename[i], fileLoci[i], mean, variance, n, numBins, threshold, upperCutoff, lowerCutoff, NSL);
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
        if(XPEHH) readAllXPEHH(filename, fileLoci, nfiles, freq1, freq2, score, chr);
        if(SOFT) readAllIHH12(filename, fileLoci, nfiles, freq1, score, chr);
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

        if(LOG_INPUT){
            cerr << "Using frequency bins from " << params.getStringFlag(ARG_LOG_INPUT) << "\n";
            flog << "Using frequency bins from " << params.getStringFlag(ARG_LOG_INPUT) << "\n";

            getMeanVarBinsFromLog(params.getStringFlag(ARG_LOG_INPUT),
                                       freq1, score, totalLoci,
                                       mean, variance,  n, numBins, threshold, (XPNSL|| XPEHH|| SOFT));
        }else{
            getMeanVarBins(freq1, score, totalLoci, mean, variance, n, numBins, threshold);
        }
        

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
            if(XPEHH) normalizeXPEHHDataByBins(filename[i], outfilename[i], fileLoci[i], mean, variance, n, numBins, threshold, upperCutoff, lowerCutoff, XPNSL);
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
int SelscanNorm::countCols(ifstream &fin)
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

int SelscanNorm::colsToSkip(ifstream &fin, int numCols)
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

void SelscanNorm::skipCols(ifstream &fin, int numCols)
{
    string junk;
    
    for(int i=0; i<numCols; i++)   
    {
        fin >> junk;
    }
}

void SelscanNorm::analyzeIHSBPWindows(string normedfiles[], int fileLoci[], int nfiles, int winSize, int numQuantiles, int minSNPs)
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
        double maxAbs = -99999.9;
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
                
                maxAbs = -99999.9;
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
            
            if(FINE_PERCENTILE){
                for(int p = 1; p <= 100; p++){
                    string perc = to_string(p) + ".0";
                    topWindowBoundary[b][perc] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 1.0 - double(p)/100.0);
                }
            }else{
                topWindowBoundary[b]["1.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.990);
                topWindowBoundary[b]["5.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.950);

            }
            
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
    
    if(FINE_PERCENTILE){
        for(int p = 1; p <= 100; p++){
            string perc = to_string(p) + ".0";
            topWindowBoundary[b][perc] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 1.0 - double(p)/100.0);
        }
    }else{
        topWindowBoundary[b]["1.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.990);
        topWindowBoundary[b]["5.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindow[start]), 1, count, 0.950);
    }

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

    //USE_GENEBED_ADDSTART
    // vector<Gene> genes;
    // if(USE_GENE_BED){
    //     readGenes(GENE_BED, genes);
    //     std::sort(genes.begin(), genes.end(), [](auto& a, auto& b){ return a.start < b.start; });
    // }
    // std::map<string, std::vector<double>> geneScores;
    // std::map<string, std::vector<double>> geneLengthsMap;
    //USE_GENEBED_ADDEND

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
        fout<<"start\tend\tn_snps\tfrac_extreme\tperc\tscore\n";  // <start> <end> <num-snps> <fraction-of-extreme-snps> <percentile> <max-score of ihs/nsl etc.>
    
        for (int j = 0; j < nSNPs[i].size(); j++)
        {
            if (nSNPs[i][j] < minSNPs || fracCrit[i][j] < 0)
            {
                fout << winStarts[i][j] << "\t" << winStarts[i][j] + winSize - 1 << "\t" << nSNPs[i][j] << "\t" << fracCrit[i][j] << "\t-1\tNA" << endl;
                continue;
            }
            double percentile = 100.0;
            for (b = 0; b < numQuantiles; b++)
            {
                if (nSNPs[i][j] <= quantileBound[b]) break;
            }

            if(FINE_PERCENTILE){
                for(int p = 1; p <= 100; p++){
                    string perc = to_string(p) + ".0";
                    if (fracCrit[i][j] >= topWindowBoundary[b][perc]){
                        percentile = double(p);
                        break;
                    }
                }
            }else{
                if (fracCrit[i][j] >= topWindowBoundary[b]["5.0"] && fracCrit[i][j] < topWindowBoundary[b]["1.0"])
                {
                    percentile = 5.0;
                }
                else if (fracCrit[i][j] >= topWindowBoundary[b]["1.0"])// && fracCrit[i][j] < topWindowBoundary[b]["0.5"])
                {
                    percentile = 1.0;
                }
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

            fout << winStarts[i][j] << "\t" << winStarts[i][j] + winSize - 1 << "\t" << nSNPs[i][j] << "\t" << fracCrit[i][j] << "\t" << percentile << "\t";
            if(maxAbsScore[i][j] == -99999.9){
                fout << "NA" << endl;
            }
            else{
                fout << maxAbsScore[i][j] << endl;
            }

            // if(!USE_GENE_BED || i != 0){
            //     fout<<"\t-"<<endl; 
            // }

            // if(USE_GENE_BED && i == 0){ // only annotate first file if multiple files   
            //     fout<<"\t";
                
            //     std::deque<Gene> active;
            //     size_t g = 0;  // pointer for genes

            //     int winStart = winStarts[i][j];
            //     int winEnd   = winStart + winSize - 1;

            //     while (!active.empty() && active.front().end < winStart) { // Remove genes that ended before this window
            //         active.pop_front();
            //     }
                
            //     while (g < genes.size() && genes[g].start <= winEnd) { // Add genes whose start is <= window end
            //         active.push_back(genes[g]);
            //         g++;
            //     }

            //     std::vector<std::string> overlaps; // Collect overlapping gene names
            //     for (const auto& gene : active) {
            //         if (gene.end >= winStart) {
            //             overlaps.push_back(gene.name);
            //             //=
            //             geneScores[gene.name].push_back(maxAbsScore[i][j]);
            //             geneLengthsMap[gene.name].push_back(gene.end - gene.start + 1);
            //             //=
            //         }
            //     }

            //     if (overlaps.empty()) {
            //         fout << "-";
            //     } else {
            //         for (size_t i = 0; i < overlaps.size(); i++) {
            //             fout << overlaps[i];
            //             if (i < overlaps.size() - 1) fout << ", ";
            //         }
            //     }
            //     fout << "\n";
            // }
        }
            
        fout.close();
        
    }

    

  
    // if(USE_GENE_BED){
    //     vector<double> geneLengths;
    //     //populate geneLengths from map
    //     for (const auto& [gene, vals] : geneLengthsMap) {
    //         if(vals.size() > 0){
    //             geneLengths.push_back(geneLengthsMap[gene][0]); // all lengths are the same for a gene
    //         }
    //     }
    //     vector<double> geneScoresMax;
    //     for (const auto& [gene, vals] : geneScores) {
    //         if(vals.size() > 0){
    //             double* data = const_cast<double*>(vals.data()); // gsl requires non-const pointer
    //             size_t n = vals.size();
    //             double maxScore  = *std::max_element(vals.begin(), vals.end());
    //             geneScoresMax.push_back(maxScore);
    //         }
    //     }
    //     geneScoresMax = regress_out_length(geneLengths, geneScoresMax);

    //     ofstream genetable; // output gene table with mean, sd, nwin, max, max_lenreg
    //     string genetablefile = GENE_BED + ".genetable";
    //     genetable.open(genetablefile.c_str());
    //     if (genetable.fail())
    //     {
    //         cerr << "ERROR: " << genetablefile << " " << strerror(errno);
    //         exit(EXIT_FAILURE);
    //     }
    //     genetable << "gene\tmean\tsd\tnwin\tmax\tmax_lenreg\n";  //gene   mean   sd   nwin   max   max_lenreg
    //     int i = 0;
    //     for (const auto& [gene, vals] : geneScores) {
    //         double* data = const_cast<double*>(vals.data()); // gsl requires non-const pointer
    //         size_t n = vals.size();
    //         double maxScore  = *std::max_element(vals.begin(), vals.end());
    //         double meanScore = gsl_stats_mean(data, 1, n);
    //         double varScore  = 0; // for n = 1
    //         if (n > 1) varScore = std::sqrt(gsl_stats_variance(data, 1, n)/n); 
    //         genetable << gene << "\t" << meanScore << "\t" << varScore << "\t" << n << "\t" << maxScore << "\t" << geneScoresMax[i] << "\t" <<"\n";
    //         i++;
    //     }
    //     genetable.close();
    // }

    delete [] quantileBound;
    delete [] topWindowBoundary;
    delete [] winStarts;
    delete [] nSNPs;
    delete [] fracCrit;
    delete [] winfilename;
    delete [] maxAbsScore; // added 

    return;
}

void SelscanNorm::analyzeXPEHHBPWindows(string normedfiles[], int fileLoci[], int nfiles, int winSize, int numQuantiles, int minSNPs)
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

            // topWindowBoundaryTop[bTop]["1.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowTop[startTop]), 1, countTop, 0.990);
            // topWindowBoundaryTop[bTop]["5.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowTop[startTop]), 1, countTop, 0.950);

            if (FINE_PERCENTILE)
            {
                for (int p = 1; p <= 100; p++)
                {
                    string perc = to_string(p) + ".0";
                    double key = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowTop[startTop]), 1, countTop, 1.0 - p / 100.0);
                    topWindowBoundaryTop[bTop][perc] = key;
                }
            }
            else
            {
                topWindowBoundaryTop[bTop]["1.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowTop[startTop]), 1, countTop, 0.990);
                topWindowBoundaryTop[bTop]["5.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowTop[startTop]), 1, countTop, 0.950);

            }



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

    if(FINE_PERCENTILE){
        for (int p = 1; p <= 100; p++)
        {
            string perc = to_string(p) + ".0";
            double key = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowTop[startTop]), 1, countTop, 1.0 - p / 100.0);
            topWindowBoundaryTop[bTop][perc] = key;
        }
    }
    else{
        topWindowBoundaryTop[bTop]["1.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowTop[startTop]), 1, countTop, 0.990);
        topWindowBoundaryTop[bTop]["5.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowTop[startTop]), 1, countTop, 0.950);
    }


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

            // topWindowBoundaryBot[bBot]["1.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowBot[startBot]), 1, countBot, 0.990);
            // topWindowBoundaryBot[bBot]["5.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowBot[startBot]), 1, countBot, 0.950);
            if(FINE_PERCENTILE){
                for(int p = 1; p<=100; p++){
                    string perc = to_string(p)+".0";
                    double key = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowBot[startBot]), 1, countBot, 1.0 - p/100.0);
                    topWindowBoundaryBot[bBot][perc] = key;
                }
            }else{
                topWindowBoundaryBot[bBot]["1.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowBot[startBot]), 1, countBot, 0.990);
                topWindowBoundaryBot[bBot]["5.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowBot[startBot]), 1, countBot, 0.950);
            }

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
    // topWindowBoundaryBot[bBot]["1.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowBot[startBot]), 1, countBot, 0.990);
    // topWindowBoundaryBot[bBot]["5.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowBot[startBot]), 1, countBot, 0.950);
    if(FINE_PERCENTILE){
        for(int p = 1; p<=100; p++){
            string perc = to_string(p)+".0";
            double key = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowBot[startBot]), 1, countBot, 1.0 - p/100.0);
            topWindowBoundaryBot[bBot][perc] = key;
        }
    }else{
        topWindowBoundaryBot[bBot]["1.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowBot[startBot]), 1, countBot, 0.990);
        topWindowBoundaryBot[bBot]["5.0"] = gsl_stats_quantile_from_sorted_data(&(allFracCritPerWindowBot[startBot]), 1, countBot, 0.950);
    }


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
        fout<<"start\tend\tnSNPs\tfrac_top\tfrac_bottom\tperc_top\tperc_bottom\ttop_score\tbottom_score\n";  
        cerr << "Creating window file " << winfilename[i] << endl; 
        flog << "Creating window file " << winfilename[i] << endl;

        for (int j = 0; j < nSNPs[i].size(); j++)
        {
            fout << winStarts[i][j] << "\t" << winStarts[i][j] + winSize - 1 << "\t" << nSNPs[i][j] << "\t" << fracCritTop[i][j] << "\t" << fracCritBot[i][j] << "\t";
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


            if(FINE_PERCENTILE){
                for(int p = 1; p <= 100; p++){
                    string perc = to_string(p) + ".0";
                    if (fracCritTop[i][j] >= topWindowBoundaryTop[bTop][perc]){
                        percentile = double(p);
                        break;
                    }
                }
            }else{

                if (fracCritTop[i][j] >= topWindowBoundaryTop[bTop]["5.0"] && fracCritTop[i][j] < topWindowBoundaryTop[bTop]["1.0"])
                {
                    percentile = 5.0;
                }
                else if (fracCritTop[i][j] >= topWindowBoundaryTop[bTop]["1.0"])// && fracCritTop[i][j] < topWindowBoundaryTop[b]["0.5"])
                {
                    percentile = 1.0;
                }
            }
            
            
            fout << percentile << "\t";

            percentile = 100.0;
            for (bBot = 0; bBot < numQuantiles; bBot++)
            {
                if (nSNPs[i][j] <= quantileBoundBot[bBot]) break;
            }

            if(FINE_PERCENTILE){
                for(int p = 1; p <= 100; p++){
                    string perc = to_string(p) + ".0";
                    if (fracCritBot[i][j] >= topWindowBoundaryBot[bBot][perc]){
                        percentile = double(p);
                        break;
                    }
                }
            }else{
                if (fracCritBot[i][j] >= topWindowBoundaryBot[bBot]["5.0"] && fracCritBot[i][j] < topWindowBoundaryBot[bBot]["1.0"])        {
                    percentile = 5.0;
                }
                else if (fracCritBot[i][j] >= topWindowBoundaryBot[bBot]["1.0"])// && fracCritTop[i][j] < topWindowBoundaryTop[b]["0.5"])
                {
                    percentile = 1.0;
                }
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

void SelscanNorm::analyzeIHH12BPWindows(string normedfiles[], int fileLoci[], int nfiles, int winSize, int numQuantiles, int minSNPs)
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

    string name, header;
    int pos;
    double gpos, freq1, ihh12, data, normedData;
    bool crit;
    int numWindows = 0;
    double maxAbs = -99999.9;
    
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
                maxAbsScore[i].push_back(maxAbs);
                maxAbs = -99999.9;

                winStart += winSize;
                winEnd += winSize;
                numSNPs = 0;
                numCrit = 0;
                
            }

            if(normedData > maxAbs) maxAbs = normedData;
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
        fout<<"start\tend\tnSNPs\tfrac_extreme\tperc\tscore\tannotation\n";  // <start> <end> <num-snps> <fraction-of-extreme-snps> <percentile> <max-score of ihs/nsl etc.>
    
        cerr << "Creating window file " << winfilename[i] << endl;
        flog << "Creating window file " << winfilename[i] << endl;
        for (int j = 0; j < nSNPs[i].size(); j++)
        {
            if (nSNPs[i][j] < minSNPs || fracCrit[i][j] < 0)
            {
                fout << winStarts[i][j] << "\t" << winStarts[i][j] + winSize - 1 << "\t" << nSNPs[i][j] << "\t" << fracCrit[i][j] << "\t-1" <<"\tNA"<< endl;
                continue;
            }
            double percentile = 100.0;
            for (b = 0; b < numQuantiles; b++)
            {
                if (nSNPs[i][j] <= quantileBound[b]) break;
            }

            if(FINE_PERCENTILE){
                for(int p = 1; p <= 100; p++){
                    string perc = to_string(p) + ".0";
                    if (fracCrit[i][j] >= topWindowBoundary[b][perc]){
                        percentile = double(p);
                        break;
                    }
                }
            }else{
                if (fracCrit[i][j] >= topWindowBoundary[b]["5.0"] && fracCrit[i][j] < topWindowBoundary[b]["1.0"])
                {
                    percentile = 5.0;
                }
                else if (fracCrit[i][j] >= topWindowBoundary[b]["1.0"])// && fracCrit[i][j] < topWindowBoundary[b]["0.5"])
                {
                    percentile = 1.0;
                }
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
            fout << winStarts[i][j] << "\t" << winStarts[i][j] + winSize - 1 << "\t" << nSNPs[i][j] << "\t" << fracCrit[i][j] << "\t" << percentile << "\t";
            if(maxAbsScore[i][j] == -99999.9){
                fout << "NA";
            }
            else{
                fout << maxAbsScore[i][j] ;
            }
            fout<<endl;
        }
        fout.close();
    }

    delete [] quantileBound;
    delete [] topWindowBoundary;
    delete [] winStarts;
    delete [] nSNPs;
    delete [] fracCrit;
    delete [] winfilename;

    delete [] maxAbsScore;

    return;
}

//Calculates mean and variance for each frequency bin
//and normalizes the data array based on these values
//freq[]: derived allele frequency array
//data[]: unnormalized score array  
void SelscanNorm::getMeanVarBins(double freq[], double data[], int nloci, double mean[], double variance[], int n[], int numBins, double threshold[])
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

void SelscanNorm::getMeanVarBinsFromLog(const std::string &binFile,
                                       double freq[], double data[], int nloci,
                                       double mean[], double variance[], int n[],
                                       int numBins, double threshold[], bool XPORSOFT)
{
    // 1. Initialize arrays
    for (int b = 0; b < numBins; b++) {
        n[b] = 0;
        mean[b] = NAN;
        variance[b] = NAN;
    }

    // 2. Open file
    std::ifstream in(binFile);
    if (!in) {
        std::cerr << "Error: cannot open bin data file " << binFile << "\n";
        exit(1);
    }

    std::string line;
    bool headerFound = false;

    // 3. Read until we find the header that starts with "bin"
    while (std::getline(in, line)) {
            if(!XPORSOFT){
                if (line.rfind("bin\tnum\mean\variance", 0) == 0) {  // starts with "bin"
                    headerFound = true;
                    break;
                }
            }else{
                if (line.rfind("num\tmean\tvariance", 0) == 0) {  // starts with "num"
                    headerFound = true;
                    break;
                }
            }
        
    }

    if (!headerFound) {
        if(XPORSOFT){
            std::cerr << "Error: no header starting with 'num\tmean\tvariance' found in " << binFile << "\n";
            exit(1);
            
        }else{
            std::cerr << "Error: no header starting with 'bin\tnum\tmean\tvariance' found in " << binFile << "\n";
            exit(1);
        }
    }

    // 4. Read the rest of the file (bin rows)
    if(!XPORSOFT){
        while (std::getline(in, line)) {
            if (line.empty()) continue;

            std::istringstream iss(line);
            double bin;
            int num;
            double m, v;
            if (!(iss >> bin >> num >> m >> v)) continue;

            // Find the right bin slot by threshold match
            for (int b = 0; b < numBins; b++) {
                if (std::abs(bin - threshold[b]) < 1e-6) {
                    n[b] = num;
                    mean[b] = m;
                    variance[b] = v;
                    break;
                }
            }
        }
    }else{
        if(numBins != 1){
            std::cerr << "Error: XP-EHH or XP-nSL or iHH12 log file should contain stats for exactly one bin.\n";
            exit(1);
        }
        int b = 0;
        while (std::getline(in, line)) {
            if (line.empty()) continue;

            std::istringstream iss(line);
            int num;
            double m, v;
            if (!(iss >> num >> m >> v)) continue;

            if (b < numBins) {
                n[b] = num;
                mean[b] = m;
                variance[b] = v;
                b++;
            }
        }
    }
    
    in.close();

    // 5. Normalize the data using the precomputed stats
    for (int i = 0; i < nloci; i++) {
        if (data[i] == MISSING) continue;
        for (int b = 0; b < numBins; b++) {
            if (freq[i] < threshold[b]) {
                if (n[b] > 1 && !std::isnan(mean[b]) && !std::isnan(variance[b])) {
                    data[i] = (data[i] - mean[b]) / std::sqrt(variance[b]);
                }
                break;
            }
        }
    }
}

//Reads a file, calculates the normalized score, and
//outputs the original row plus normed score
void SelscanNorm::normalizeIHSDataByBins(string &filename, string &outfilename, int &fileLoci, double mean[], double variance[], int n[], int numBins, double threshold[], double upperCutoff, double lowerCutoff, bool NSL)
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
    string chr;
    double freq, data, normedData, ihh1, ihh2;;
    int numInBin = 0;
    string junk;
    getline(fin, header);

    const std::string expectedHeaderIHS = "chr\tid\tpos\tfreq\tihh1\tihh0\tihs";
    const std::string expectedHeaderNSL = "chr\tid\tpos\tfreq\tihh1\tihh0\tnsl";


    if (NSL){
         if(expectedHeaderNSL != header){
            cerr << "ERROR: Expected header '" << expectedHeaderNSL << "' but found '" << header << "' in file " << filename << "\n";
            flog << "ERROR: Expected header '" << expectedHeaderNSL << "' but found '" << header << "' in file " << filename << "\n";
            exit(EXIT_FAILURE);
        }
        fout << header + "\tnorm_nsl\tcrit\n"; 
    } else{
         if(expectedHeaderIHS != header){
            cerr << "ERROR: Expected header '" << expectedHeaderIHS <<  "' but found '" << header << "' in file " << filename << "\n";
            flog << "ERROR: Expected header '" << expectedHeaderIHS <<  "' but found '" << header << "' in file " << filename << "\n";
            exit(EXIT_FAILURE);
        }
        fout << header + "\tnorm_ihs\tcrit\n";  
    } 

    // determine the number of cols to skip 
    // (the number more than the number we care about: 6)
    int numColsToSkip = 0;
    numColsToSkip = colsToSkip(fin, 6);

    for (int j = 0; j < fileLoci; j++)
    {
        fin >> chr;
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
            fout << chr << "\t" << name << "\t"
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

void SelscanNorm::normalizeXPEHHDataByBins(string &filename, string &outfilename, int &fileLoci, double mean[], double variance[], int n[], int numBins, double threshold[], double upperCutoff, double lowerCutoff, bool XPNSL)
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
    string chr;
    double gpos, freq1, freq2, data, normedData, ihh1, ihh2;;
    int numInBin = 0;

    getline(fin, header);
    if (XPNSL){
        fout << header + "\tnormxpnsl\tcrit\n"; 
    } else{
        fout << header + "\tnormxpehh\tcrit\n";  
    } 

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


void SelscanNorm::normalizeIHH12DataByBins(string &filename, string &outfilename, int &fileLoci, double mean[], double variance[], int n[], int numBins, double threshold[], double upperCutoff, double lowerCutoff)
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

    const std::string expectedHeader = "chr\tid\tpos\tp1\tihh12";
    if (header != expectedHeader) {     //check that header matches expected format
        flog << "ERROR: Expected header '" << expectedHeader << "' but found '" << header << "' in file " << filename << "\n";
        cerr << "ERROR: Expected header '" << expectedHeader << "' but found '" << header << "' in file " << filename << "\n";
        exit(EXIT_FAILURE);
    }

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
int SelscanNorm::checkIHSfile(ifstream &fin)
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

    nloci--; //added in v3, skip header line

    fin.clear();
    fin.seekg(start);

    return nloci;
}

void SelscanNorm::readAllIHS(vector<string> filename, int fileLoci[], int nfiles, double freq[], double score[], string& chr)
{
    ifstream fin;
    string junk;
    int overallCount = 0;
    for (int i = 0; i < nfiles; i++)
    {
        fin.open(filename[i].c_str());
        getline(fin, junk); //added v3, skip header line

        int numColsToSkip = 0;
        numColsToSkip = colsToSkip(fin, 6);

        for (int j = 0; j < fileLoci[i]; j++)
        {
            fin >> chr;
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

int SelscanNorm::checkXPEHHfile(ifstream &fin)
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


int SelscanNorm::checkIHH12file(ifstream &fin)
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

void SelscanNorm::readAllXPEHH(vector<string> filename, int fileLoci[], int nfiles, double freq1[], double freq2[], double score[], string& chr)
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
            fin >> chr;
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

void SelscanNorm::readAllIHH12(vector<string> filename, int fileLoci[], int nfiles, double freq1[], double score[], string& chr)
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
            fin >> chr;
            fin >> junk;
            fin >> freq1[overallCount];
            fin >> score[overallCount];
            overallCount++;
        }
        fin.close();
    }

    return;
}


int SelscanNorm::countFields(const string &str)
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

bool SelscanNorm::isint(string str)
{
    for (string::iterator it = str.begin(); it != str.end(); it++)
    {
        if (!isdigit(*it)) return 0;
    }

    return 1;
}


std::pair<double,double> SelscanNorm::fit_length_regression(
    const std::vector<double> &lengths,
    const std::vector<double> &scores)
{
    size_t n = lengths.size();
    size_t p = 2; // intercept + log(length)

    // Collect valid entries
    std::vector<size_t> valid_idx;
    valid_idx.reserve(n);
    for (size_t i = 0; i < n; i++) {
        if (abs(scores[i]) != MISSING_SCORE)
            valid_idx.push_back(i);
    }

    size_t m = valid_idx.size();
    if (m == 0) return {0.0, 0.0}; // nothing to fit

    gsl_matrix *X = gsl_matrix_alloc(m, p);
    gsl_vector *y = gsl_vector_alloc(m);
    gsl_vector *c = gsl_vector_alloc(p);
    gsl_matrix *cov = gsl_matrix_alloc(p, p);
    double chisq;

    for (size_t j = 0; j < m; j++) {
        size_t i = valid_idx[j];
        gsl_matrix_set(X, j, 0, 1.0);
        gsl_matrix_set(X, j, 1, std::log(lengths[i]));
        gsl_vector_set(y, j, scores[i]);
    }

    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(m, p);
    gsl_multifit_linear(X, y, c, cov, &chisq, work);
    gsl_multifit_linear_free(work);

    double beta0 = gsl_vector_get(c, 0);
    double beta1 = gsl_vector_get(c, 1);

    gsl_matrix_free(X);
    gsl_matrix_free(cov);
    gsl_vector_free(y);
    gsl_vector_free(c);

    return {beta0, beta1};
}

std::vector<double> SelscanNorm::compute_length_residuals(
    const std::vector<double> &lengths,
    const std::vector<double> &scores,
    double beta0,
    double beta1)
{
    size_t n = lengths.size();
    std::vector<double> residuals(n, -1.0);

    for (size_t i = 0; i < n; i++) {
        if (abs(scores[i]) != MISSING_SCORE) {
            double pred = beta0 + beta1 * std::log(lengths[i]);
            residuals[i] = scores[i] - pred;
        }
    }

    return residuals;
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
std::vector<double> SelscanNorm::regress_out_length(
    const std::vector<double> &lengths,
    const std::vector<double> &scores) 
{
    size_t n = lengths.size();
    size_t p = 2; // intercept + log(length)

    // First collect only valid (score != -1) entries
    std::vector<size_t> valid_idx;
    valid_idx.reserve(n);
    for (size_t i = 0; i < n; i++) {
        if (abs(scores[i]) != MISSING_SCORE) {
            valid_idx.push_back(i);
        }
    }

    size_t m = valid_idx.size(); // number of valid points
    if (m == 0) {
        // nothing to fit, just return original
        return scores;
    }

    gsl_matrix *X = gsl_matrix_alloc(m, p);
    gsl_vector *y = gsl_vector_alloc(m);
    gsl_vector *c = gsl_vector_alloc(p);
    gsl_matrix *cov = gsl_matrix_alloc(p, p);
    double chisq;

    // Fill design matrix and response with only valid data
    for (size_t j = 0; j < m; j++) {
        size_t i = valid_idx[j];
        gsl_matrix_set(X, j, 0, 1.0);                  // intercept
        gsl_matrix_set(X, j, 1, std::log(lengths[i])); // log(length)
        gsl_vector_set(y, j, scores[i]);
    }

    // Fit linear regression
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(m, p);
    gsl_multifit_linear(X, y, c, cov, &chisq, work);
    gsl_multifit_linear_free(work);

    double beta0 = gsl_vector_get(c, 0);
    double beta1 = gsl_vector_get(c, 1);

    // Compute residuals (keep same size as input)
    std::vector<double> residuals(n, -1.0);
    for (size_t j = 0; j < m; j++) {
        size_t i = valid_idx[j];
        double pred = beta0 + beta1 * std::log(lengths[i]);
        residuals[i] = scores[i] - pred;
    }

    // Shift residuals so the minimum valid one is 0
    // double minResidual = std::numeric_limits<double>::max();
    // for (double r : residuals) {
    //     if (r != -1.0) minResidual = std::min(minResidual, r);
    // }
    // if (minResidual < 0.0) {
    //     for (auto &r : residuals) {
    //         if (r != -1.0) r -= minResidual;
    //     }
    // }

    // Free memory
    gsl_matrix_free(X);
    gsl_matrix_free(cov);
    gsl_vector_free(y);
    gsl_vector_free(c);

    return residuals;
}

// std::vector<double> SelscanNorm::regress_out_length(const std::vector<double> &lengths,
//                                        const std::vector<double> &scores) {
//     size_t n = lengths.size();
//     size_t p = 2; // intercept + log(length)

    
//     gsl_matrix *X = gsl_matrix_alloc(n, p);
//     gsl_vector *y = gsl_vector_alloc(n);
//     gsl_vector *c = gsl_vector_alloc(p);
//     gsl_matrix *cov = gsl_matrix_alloc(p, p);
//     double chisq;

//     // Fill design matrix and response
//     for (size_t i = 0; i < n; i++) {
//         gsl_matrix_set(X, i, 0, 1.0);                  // intercept
//         gsl_matrix_set(X, i, 1, std::log(lengths[i])); // log(length)
//         gsl_vector_set(y, i, scores[i]);
//     }

//     // Fit linear regression
//     gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n, p);
//     gsl_multifit_linear(X, y, c, cov, &chisq, work);
//     gsl_multifit_linear_free(work);

//     double beta0 = gsl_vector_get(c, 0);
//     double beta1 = gsl_vector_get(c, 1);

//     // Compute residuals
//     std::vector<double> residuals(n);
//     for (size_t i = 0; i < n; i++) {
//         double pred = beta0 + beta1 * std::log(lengths[i]);
//         residuals[i] = scores[i] - pred;
//     }

//     //After computing residuals, shift them so the minimum is 0
//     double minResidual = *std::min_element(residuals.begin(), residuals.end());
//     if (minResidual < 0.0) {
//         for (auto &r : residuals) r -= minResidual;
//     }
//     // Free memory
//     gsl_matrix_free(X);
//     gsl_matrix_free(cov);
//     gsl_vector_free(y);
//     gsl_vector_free(c);

    
//     return residuals;
// }