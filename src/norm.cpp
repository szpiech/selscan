#include <algorithm>
#include "norm.h"
#include "normio.h"
#include <limits>
#include <iomanip>
#include <deque>

// #include <gsl/gsl_statistics_double.h>
// #include <gsl/gsl_interp.h>

//TODO
// TODO: winfile allow multiple windows

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
    params.addFlag(ARG_WINOFFSET, DEFAULT_WINOFFSET, "SILENT", HELP_WINOFFSET);



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
    params.addFlag(ARG_GTF, DEFAULT_GTF, "", HELP_GTF); //added in v3.1
    params.addFlag(ARG_GENE_OUT, DEFAULT_GENE_OUT, "", HELP_GENE_OUT); // U15

    params.addListFlag(ARG_WIN_FILES, DEFAULT_WIN_FILES, "", HELP_WIN_FILES);
    params.addListFlag(ARG_NORM_FILES, DEFAULT_NORM_FILES, "", HELP_NORM_FILES); //added in v3
    
    //params.addFlag(ARG_WIN_FILE, DEFAULT_WIN_FILE, "", HELP_WIN_FILE); //added in v3
    
    params.addFlag(ARG_GENE_SETA, DEFAULT_GENE_SETA, "", HELP_GENE_SETA); //added in v3
    params.addFlag(ARG_GENE_SETB, DEFAULT_GENE_SETB, "", HELP_GENE_SETB); //added in v3
    params.addFlag(ARG_PERM_SEED, DEFAULT_PERM_SEED, "", HELP_PERM_SEED);

    params.addFlag(ARG_LOG_ONLY, DEFAULT_LOG_ONLY, "", HELP_LOG_ONLY); //added in v3
    // params.addFlag(ARG_NO_HEADER, DEFAULT_NO_HEADER, "", HELP_NO_HEADER); //added in v3

    try
    {
        params.parseCommandLine(argc, argv);
    }
    catch (...)
    {
        return 1;
    }

    this->WINOFFSET = params.getIntFlag(ARG_WINOFFSET);
    int numBins = params.getIntFlag(ARG_FREQ_BINS);
    vector<string> filename = params.getStringListFlag(ARG_FILES);
    
    const int winSize = params.getIntFlag(ARG_WINSIZE);
    string infoOutfile = params.getStringFlag(ARG_LOG);
    const int numQBins = params.getIntFlag(ARG_QBINS);
    const int minSNPs = params.getIntFlag(ARG_MINSNPS);
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

    // read window file
    vector<string> windowFiles = params.getStringListFlag(ARG_WIN_FILES);
    vector<string> normFiles = params.getStringListFlag(ARG_NORM_FILES);

    //string windowFile = params.getStringFlag(ARG_WIN_FILE);
    
    string geneBedFile = params.getStringFlag(ARG_BED);
    string geneGTFFile = params.getStringFlag(ARG_GTF);
    string geneOut = params.getStringFlag(ARG_GENE_OUT);

    //bool HAVE_WINFILE; //= params.flagWasSet(ARG_WIN_FILES);
    //bool HAVE_NORMFILE; // = params.flagWasSet(ARG_NORM_FILES);

    string geneSetA = params.getStringFlag(ARG_GENE_SETA);
    string geneSetB = params.getStringFlag(ARG_GENE_SETB);

    
    bool PERMUTE_TEST = (geneSetA != DEFAULT_GENE_SETA) && (geneSetB != DEFAULT_GENE_SETB);
    

    GeneAnalyzer genex;
    //--win-files + --gene-bed 
    //--files + --gene-bed/gtf/gff3
    //--norm-files + --gene-bed/gtf/gff3 + --bp-win


    bool useGTF = false;
    bool useBED = false;
    if(geneGTFFile != DEFAULT_GTF) useGTF = true;
    if(geneBedFile != DEFAULT_BED) useBED = true;
    
    //check that both gtf and bed arent provided
    if(useGTF && useBED){
        cerr << "ERROR: Cannot provide both " + ARG_BED + " and " + ARG_GTF + " Use one for gene annotations.\n";
        return 1;
    }

    string geneFile = useGTF ? geneGTFFile : geneBedFile;

    bool LOG_INPUT =  params.flagWasSet(ARG_LOG_INPUT);
    bool LOG_ONLY = params.getBoolFlag(ARG_LOG_ONLY);
    
    if(params.getStringFlag(ARG_LOG_INPUT) == infoOutfile){
        cerr << "ERROR: --log-input and --log cannot have the same filename: both are " << infoOutfile << "\n";
        exit(EXIT_FAILURE);
    }
    
    if(PERMUTE_TEST){
        cerr << "Performing permutation test between gene sets.\n";
        // if(geneSetA == DEFAULT_GENE_SETA || geneSetB == DEFAULT_GENE_SETB){
        //     cerr << "ERROR: --gene-target and --gene-background must be provided to run permutation test.\n";
        //     return 1;
        // }
        cerr << "Running permutation test with target genes from " << geneSetA << " and background genes from " << geneSetB << "\n";
        genex.perm_test(geneSetA, geneSetB, params.getIntFlag(ARG_PERM_SEED));
        return 0;
    }

    // B20 (PI decision Q11): with --log-input the number of frequency bins is
    // the number of rows in the log's bin table (3.0.1 always used 100, so
    // logs written with other --bins silently dropped most loci). --bins may
    // be given only if it agrees.
    if (LOG_INPUT && (IHS || NSL)) {
        int nLogBins = countLogBins(params.getStringFlag(ARG_LOG_INPUT));
        if (nLogBins <= 0) {
            cerr << "ERROR: no frequency-bin table (header 'bin\tnum\tmean\tvariance') in "
                 << params.getStringFlag(ARG_LOG_INPUT) << ".\n";
            exit(1);
        }
        if (params.flagWasSet(ARG_FREQ_BINS) && numBins != nLogBins) {
            cerr << "ERROR: --bins " << numBins << " disagrees with the " << nLogBins
                 << " frequency bins in " << params.getStringFlag(ARG_LOG_INPUT) << ".\n";
            exit(1);
        }
        numBins = nLogBins;
    }

    if(params.flagWasSet(ARG_FREQ_BINS) && (params.flagWasSet(ARG_XPEHH) || params.flagWasSet(ARG_XPNSL) || params.flagWasSet(ARG_SOFT))){
        cerr << "ERROR: Options --bins and --xpehh/--xpnsl/--ihh12 cannot be used together. "
         << "Frequency bins are only used for iHS and nSL normalization.\n";
        exit(1);
    }
       
    if(params.flagWasSet(ARG_QBINS) && LOG_INPUT){
        cerr << "ERROR: Options --log-input and --qbins cannot be used together. "
         << "--log-input already provides quantile bin information.\n";
        exit(1);
    }

    

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

    // B35: with 0, empty windows were counted but never filled, and their
    // uninitialised values went into the quantile tables
    if (minSNPs < 1)
    {
        cerr << "ERROR: " << ARG_MINSNPS << " must be at least 1.\n";
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


    if(XPNSL) XPEHH = true;


    // cout<<("=== Parameters ===\n");
    // cout<<params.flagWasSet(ARG_FILES)<<endl;
    
    if( filename[0] == DEFAULT_FILES &&  normFiles[0]==DEFAULT_NORM_FILES &&  windowFiles[0]==DEFAULT_WIN_FILES ){
        cerr << "ERROR: Must provide at least one input file with " + ARG_FILES + " for normalization.\n";
        return 1;
    }


    bool HAVE_FILES = (filename[0]!=DEFAULT_FILES); // user gave files (e.g., .ihs.out, .nsl.out, etc.)
    bool HAVE_NORMFILE = (normFiles[0]!=DEFAULT_NORM_FILES); // user gave normfile
    bool HAVE_WINFILE = (windowFiles[0]!=DEFAULT_WIN_FILES); // user gave winfile

    bool DO_NORM = HAVE_FILES;                       // if input stat files (.out) provided
    bool DO_GENE = useBED || useGTF;                // user gave --gene-bed (.bed)

    int nfiles;
    if(HAVE_WINFILE){
        nfiles = windowFiles.size();
    }else if(HAVE_NORMFILE){
        nfiles = normFiles.size();
    }else{
        nfiles = filename.size();
    }

    //cant prode normfile winfile with files
    if( (HAVE_NORMFILE || HAVE_WINFILE) && HAVE_FILES ){
        cerr << "ERROR: Cannot provide input files with " + ARG_FILES + " when providing " + ARG_NORM_FILES + " or " + ARG_WIN_FILES + ".\n";
        return 1;
    }

    //cant provide winfile with normfile
    if( HAVE_NORMFILE && HAVE_WINFILE ){
        cerr << "ERROR: Cannot provide both " + ARG_NORM_FILES + " and " + ARG_WIN_FILES + ".\n";
        return 1;
    }

    // if --win-files provided, --bp-win is not meaningful
    if( BPWIN && HAVE_WINFILE ){
        cerr << "ERROR: Cannot use " + ARG_BPWIN + " with " + ARG_WIN_FILES + ".\n";
        return 1;
    }

    if( HAVE_NORMFILE ){
        nfiles = normFiles.size();
        filename = normFiles;
    }

    if( HAVE_WINFILE ){
        nfiles = windowFiles.size();
        filename = windowFiles;
    }


    // Allocated on the --files (DO_NORM) or --norm-files --bp-win paths only;
    // freed once at the end (B23).
    string *outfilename = nullptr;
    int *fileLoci = nullptr;
    
    // B49: the log is opened for every mode (3.0.1 opened it only with
    // --files and closed it before the window and gene analyses, so their
    // tables and warnings never reached it; --norm-files wrote no log).
    flog.open(infoOutfile.c_str());
    if (flog.fail())
    {
        cerr << "ERROR: Cannot open " << infoOutfile << " for logging " << strerror(errno) << endl;
        return 1;
    }
    for (int i = 0; i < argc; i++)
    {
        flog << argv[i] << " ";
    }
    flog << "\n\n";

    // column name of the score, used to check input headers (U3); --xpnsl also sets XPEHH (U14)
    const string statLabel = IHS ? "ihs" : (NSL ? "nsl" : (XPNSL ? "xpnsl" : (XPEHH ? "xpehh" : "ihh12")));

    if(DO_NORM){
        cerr << "You have provided " << nfiles << " output files for joint normalization.\n";

        outfilename = new string[nfiles];
        fileLoci = new int[nfiles];
        ifstream fin;

        int totalLoci = 0;



        //flog << "Input files:\n";

        //For each file, open it, and check it for integrity
        //Also record total number of lines so we can allocate
        //enough space for the array of paired data that will
        //be used to calculate E[X] and E[X^2]
        for (int i = 0; i < nfiles; i++)
        {
            char str[16];
            snprintf(str, sizeof str, "%d", numBins);
            if(IHS || NSL) outfilename[i] = filename[i] + "." + str + "bins.norm";
            if(XPEHH || SOFT) outfilename[i] = filename[i] + ".norm";
            
            fin.open(filename[i].c_str());
            if (fin.fail())
            {
                cerr << "ERROR: " << "Cannot open " << filename[i] << " " << endl; // << strerror(errno);
                flog << "ERROR: " << "Cannot open " << filename[i] << " " << endl; // << strerror(errno);
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
                fileLoci[i] = checkStatFile(fin, filename[i], statLabel, false);
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
        if (totalLoci == 0) {   // U11
            cerr << "ERROR: the input files have no data rows.\n";
            flog << "ERROR: the input files have no data rows.\n";
            return 1;
        }


        string chr;
        //@NORMLOGINPUT
        if (IHS || NSL)
        {
            cerr << "Reading all data.\n";
            double *freq = new double[totalLoci];
            double *score = new double[totalLoci];
            //read in all data
            readAllIHS(filename, fileLoci, nfiles, freq, score);
            reportBadRows();

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
                int nFreqOne = getMeanVarBins(freq, score, totalLoci, mean, variance, n, numBins, threshold, false);
                if (nFreqOne > 0) {
                    cerr << "Excluded " << nFreqOne << " loci with derived allele frequency 1.0 from normalisation (no frequency bin); they are not written to the .norm output.\n";
                    flog << "Excluded " << nFreqOne << " loci with derived allele frequency 1.0 from normalisation (no frequency bin); they are not written to the .norm output.\n";
                }

            }

            double upperCutoff, lowerCutoff;
            critCutoffs(score, totalLoci, critPercent, critNum, false, upperCutoff, lowerCutoff);
            delete [] freq;
            delete [] score;

            //Output bins info to file.
            cerr << "bin\tnum\tmean\tvariance\n";
            flog << "bin\tnum\tmean\tvariance\n";
            for (int i = 0; i < numBins; i++)
            {
                cerr << threshold[i] << "\t" << n[i] <<  "\t" << mean[i] << "\t" << variance[i] << endl;
                // full precision so --log-input reproduces this run exactly (B20)
                flog << threshold[i] << "\t" << n[i] <<  "\t" << std::setprecision(17) << mean[i] << "\t" << variance[i] << std::setprecision(6) << endl;
            }


            //Read each file and create normed files (not with --log-only, U2).
            if(!LOG_ONLY){
                if (FIRST){
                    cerr << "Only normalizing first file as --first flag is set.\n";
                    nfiles = 1;
                }
                for (int i = 0; i < nfiles; i++)
                {
                    cerr << "Normalizing " << filename[i] << "\n";
                    normalizeIHSDataByBins(filename[i], outfilename[i], fileLoci[i], mean, variance, n, numBins, threshold, upperCutoff, lowerCutoff, NSL);
                }
            }else{
                cerr << "Skipping outputting normalized output as --log-only flag is set.\n";
            }



            
            
        }
        else if (XPEHH || SOFT) {

            cerr << "Reading all data.\n";
            double *freq1 = new double[totalLoci];
            double *freq2 = nullptr;
            if(XPEHH) freq2 = new double[totalLoci];
            double *score = new double[totalLoci];
            //read in all data
            if(XPEHH) readAllXPEHH(filename, fileLoci, nfiles, freq1, freq2, score);
            if(SOFT) readAllIHH12(filename, fileLoci, nfiles, freq1, score);
            reportBadRows();
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
                getMeanVarBins(freq1, score, totalLoci, mean, variance, n, numBins, threshold, true);
            }
            

            double upperCutoff, lowerCutoff;
            critCutoffs(score, totalLoci, critPercent, critNum, SOFT, upperCutoff, lowerCutoff);
            delete [] freq1;
            if (XPEHH) delete [] freq2;
            delete [] score;

            //Output bins info to file.
            cerr << "num\tmean\tvariance\n";
            flog << "num\tmean\tvariance\n";
            for (int i = 0; i < numBins; i++)
            {
                cerr << n[i] <<  "\t" << mean[i] << "\t" << variance[i] << endl;
                flog << n[i] <<  "\t" << std::setprecision(17) << mean[i] << "\t" << variance[i] << std::setprecision(6) << endl;
            }

            // Only normalizing if necessary
            if(!LOG_ONLY){
                
                if (FIRST){
                    cerr << "Only normalizing first file as --first flag is set.\n";
                    nfiles = 1;
                }

                //Read each file and create normed files.
                for (int i = 0; i < nfiles; i++)
                {
                    cerr << "Normalizing " << filename[i] << "\n";
                    if(XPEHH) normalizeXPEHHDataByBins(filename[i], outfilename[i], fileLoci[i], mean, variance, n, numBins, threshold, upperCutoff, lowerCutoff, XPNSL);
                    if(SOFT) normalizeIHH12DataByBins(filename[i], outfilename[i], fileLoci[i], mean, variance, n, numBins, threshold, upperCutoff, lowerCutoff);
                    //fin[i].close();
                    //fout[i].close();
                }
            }else{
                cerr << "Skipping outputting normalized output as --log-only flag is set.\n";
            }
            
            delete [] threshold;
            delete [] mean;
            delete [] variance;
            delete [] n;

           
        }

        // delete [] threshold;
        // delete [] mean;
        // delete [] variance;
        // delete [] n;
    }




    // --log-only writes no .norm files, so there is nothing for the window
    // or gene analyses to read (U2)
    if (LOG_ONLY && DO_NORM && (BPWIN || DO_GENE)) {
        cerr << "Skipping window and gene analyses as --log-only flag is set.\n";
        flog << "Skipping window and gene analyses as --log-only flag is set.\n";
        BPWIN = false;
        DO_GENE = false;
    }

    if(BPWIN){ // start from .ihs.out, then require .norm files
        if(DO_NORM){
            cerr << "Using " << nfiles << " normalized output files for window-based analysis.\n";
        }else{
            cerr << "You have provided " << nfiles << " normalized files for window-based analysis.\n";
            outfilename = new string[nfiles];
            fileLoci = new int[nfiles];
            ifstream fin;
            int totalLoci = 0;
                   //For each file, open it, and check it for integrity
            //Also record total number of lines so we can allocate
            //enough space for the array of paired data that will
            //be used to calculate E[X] and E[X^2]
            for (int i = 0; i < nfiles; i++)
            {
                fin.open(filename[i].c_str());
                if (fin.fail())
                {
                    cerr << "ERROR: " << filename[i] << " " << strerror(errno);
                    flog << "ERROR: " << filename[i] << " " << strerror(errno);
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
                    fileLoci[i] = checkStatFile(fin, filename[i], statLabel, true);
                    totalLoci += fileLoci[i];
                }
                catch (...)
                {
                    return 1;
                }
                fin.close();
                outfilename[i] = filename[i]; //already .norm files
            }

            cerr << "\nTotal loci: " << totalLoci << endl; // across all files
            flog << "\nTotal loci: " << totalLoci << endl; // across all files
            if (totalLoci == 0) {   // U11
                cerr << "ERROR: the input files have no data rows.\n";
                flog << "ERROR: the input files have no data rows.\n";
                return 1;
            }
        }


        if(IHS || NSL){
            analyzeIHSBPWindows(outfilename, fileLoci, nfiles, winSize, numQBins, minSNPs); //require .norm files
        }
        if (XPEHH) analyzeXPEHHBPWindows(outfilename, fileLoci, nfiles, winSize, numQBins, minSNPs);
        if (SOFT) analyzeIHH12BPWindows(outfilename, fileLoci, nfiles, winSize, numQBins, minSNPs);
    }


    bool XP = (XPEHH || XPNSL);

    if(LOG_ONLY) DO_GENE = false;
    if(geneOut != DEFAULT_GENE_OUT && (!DO_GENE || HAVE_WINFILE)){
        cerr << "WARNING: " << ARG_GENE_OUT << " is ignored: no gene table is written in this run.\n";
        flog << "WARNING: " << ARG_GENE_OUT << " is ignored: no gene table is written in this run.\n";
    }
    if(DO_GENE){ 
        cerr << "Annotating " << nfiles << " files using " << geneFile << "\n";
        
        vector<string> normFilenames(nfiles);
        vector<string> winFilenames(nfiles);

        if(HAVE_WINFILE){ //--win-files
            cerr << "You have provided " << nfiles << " window files for gene annotation.\n";
            genex.annotateWindows(geneFile, useGTF, filename);
        }else if(DO_NORM){ // --files but not --win-files
            // nfiles, not filename.size(): --first sets nfiles = 1 and only
            // the first file is normalised (B24)
            cerr << "Using " << nfiles << " normalized output files for gene-based analysis.\n";
            for (int i = 0; i < nfiles; i++)
            {
                normFilenames[i] = outfilename[i];
            }
        }else{ // no --files, no --win-files, but maybe --norm-files
            cerr << "You have provided " << nfiles << " normalized files for gene-based analysis.\n";
            for (int i = 0; i < nfiles; i++)
            {
                // one stream per file: reopening an open ifstream fails, so
                // 3.0.1 checked only the first file (U3)
                ifstream fin(filename[i].c_str());
                if (fin.fail()) {
                    cerr << "ERROR: Cannot open " << filename[i] << "\n";
                    flog << "ERROR: Cannot open " << filename[i] << "\n";
                    return 1;
                }
                try { checkStatFile(fin, filename[i], statLabel, true); }
                catch (...) { return 1; }
                normFilenames[i] = normFiles[i];
            }
            
        }

        if(BPWIN) { // must be with --norm-files or --files
            // the window tables written by analyze*BPWindows, not the .norm
            // files they were built from (B48)
            for(int i = 0; i < nfiles; i++)
            {
                winFilenames[i] = outfilename[i] + "." + to_string(winSize / 1000) + "kb.windows";
            }
            genex.annotateWindows(geneFile, useGTF, winFilenames);
            
        }


        if(!HAVE_WINFILE){
            string stat = IHS ? "ihs" : (NSL ? "nsl" : (XPNSL ? "xpnsl" : (XPEHH ? "xpehh" : "ihh12")));   // --xpnsl also sets XPEHH (U14)
            // U15: the gene table goes next to the first normalized input,
            // named after it and the gene file, so runs that share one BED
            // do not overwrite each other (3.0.1 wrote <bed>.<stat>.genetable
            // beside the BED). --gene-out sets the path.
            string tableFile = geneOut;
            if (tableFile == DEFAULT_GENE_OUT) {
                size_t slash = geneFile.find_last_of("/\\");
                string geneBase = (slash == string::npos) ? geneFile : geneFile.substr(slash + 1);
                tableFile = normFilenames[0] + "." + geneBase + ".genetable";
            }
            genex.annotateSNPs(geneFile, useGTF, normFilenames, XP, stat, minSNPs, tableFile);
            
        }
        
    }
            
    delete [] outfilename;
    delete [] fileLoci;
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

  char str[16];
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

// Window percentile tables for --bp-win (B22, U4).
// Windows (sorted by SNP count) are split into numQuantiles SNP-count bins
// with the bounds quantileBound[b] = the (b+1)/numQuantiles quantile of the
// SNP counts. Each window goes to the first bin with nSNPs <= bound, the
// same rule the output uses (3.0.1 built the slices with a different rule:
// slices were shifted by one window and, with tied bounds, some bins were
// never filled, so their windows were compared with 0 and labelled top 1%).
// Within a bin, the 99th and 95th percentiles of the fraction of extreme
// SNPs (or every percentile with --fine-percentile) are recorded.
void SelscanNorm::windowQuantiles(const double sortedSNPs[], double frac[], int numWindows, int numQuantiles,
                                  double quantileBound[], map<string, double> topWindowBoundary[], const string &label)
{
    for (int q = 0; q < numQuantiles; q++)
        quantileBound[q] = (numWindows > 0)
            ? snum::quantile_sorted(sortedSNPs, numWindows, double(q + 1) / double(numQuantiles))
            : 0.0;
    if (numWindows == 0) {
        cerr << "WARNING: no windows" << label << " with enough SNPs; no window percentiles computed.\n";
        flog << "WARNING: no windows" << label << " with enough SNPs; no window percentiles computed.\n";
        return;
    }

    vector<int> qStart(numQuantiles, 0), qCount(numQuantiles, 0);
    int bb = 0;
    for (int i = 0; i < numWindows; i++) {
        while (bb < numQuantiles - 1 && sortedSNPs[i] > quantileBound[bb]) bb++;
        if (qCount[bb] == 0) qStart[bb] = i;
        qCount[bb]++;
    }

    cerr << "\n" << label << (label.empty() ? "" : "\n") << "nSNPs 1.0 5.0 nWindows\n";
    flog << "\n" << label << (label.empty() ? "" : "\n") << "nSNPs 1.0 5.0 nWindows\n";
    int nSmall = 0, nUsed = 0;
    for (bb = 0; bb < numQuantiles; bb++) {
        if (qCount[bb] == 0) continue;   // tied bounds: no window is assigned here
        nUsed++;
        double *slice = &(frac[qStart[bb]]);
        snum::sort(slice, qCount[bb]);
        if (FINE_PERCENTILE) {
            for (int p = 1; p <= 100; p++)
                topWindowBoundary[bb][to_string(p) + ".0"] = snum::quantile_sorted(slice, qCount[bb], 1.0 - double(p) / 100.0);
        } else {
            topWindowBoundary[bb]["1.0"] = snum::quantile_sorted(slice, qCount[bb], 0.990);
            topWindowBoundary[bb]["5.0"] = snum::quantile_sorted(slice, qCount[bb], 0.950);
        }
        if (qCount[bb] < 100) nSmall++;
        cerr << quantileBound[bb] << " " << topWindowBoundary[bb]["1.0"] << " " << topWindowBoundary[bb]["5.0"] << " " << qCount[bb] << "\n";
        flog << quantileBound[bb] << " " << topWindowBoundary[bb]["1.0"] << " " << topWindowBoundary[bb]["5.0"] << " " << qCount[bb] << "\n";
    }
    cerr << "\n"; flog << "\n";
    if (nSmall > 0) {
        cerr << "WARNING: " << nSmall << " of " << nUsed << " SNP-count bins" << label << " hold fewer than 100 windows; "
             << "their 1% and 5% window labels rest on few windows (use fewer --qbins or more data).\n";
        flog << "WARNING: " << nSmall << " of " << nUsed << " SNP-count bins" << label << " hold fewer than 100 windows; "
             << "their 1% and 5% window labels rest on few windows (use fewer --qbins or more data).\n";
    }
}

// Percentile label of a window: the first percentile whose within-bin
// boundary the window's fraction of extreme SNPs EXCEEDS (U4, PI decision
// Q10: '>' rather than '>='; with '>=' every window tied at the boundary,
// often 0, was labelled 1%). 100 means not in the top 5% (or any fine
// percentile).
double SelscanNorm::windowPercentile(int nSNPs, double frac, const double quantileBound[], int numQuantiles,
                                     map<string, double> topWindowBoundary[])
{
    int b = 0;
    while (b < numQuantiles - 1 && nSNPs > quantileBound[b]) b++;
    const map<string, double> &m = topWindowBoundary[b];
    if (m.empty()) return 100.0;
    if (FINE_PERCENTILE) {
        for (int p = 1; p <= 100; p++)
            if (frac > m.at(to_string(p) + ".0")) return double(p);
        return 100.0;
    }
    if (frac > m.at("1.0")) return 1.0;
    if (frac > m.at("5.0")) return 5.0;
    return 100.0;
}

void SelscanNorm::analyzeIHSBPWindows(string normedfiles[], int fileLoci[], int nfiles, int winSize, int numQuantiles, int minSNPs)
{
    vector<string> chrs(nfiles);
    
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

    char str[16];
    snprintf(str, sizeof str, "%d", winSize / 1000);

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
        int winStart = 1 + this->WINOFFSET;
        int winEnd = winStart + winSize - 1;
        // B36, B46: the window analysis needs one chromosome per file in
        // position order, and starts at --win-offset
        const int firstWinStart = winStart;
        string prevChr;
        int prevPos = 0;
        int numSNPs = 0;
        int numCrit = 0;
        double maxAbs = -99999.9;

        //skip header line
        string header, junk;
        getline(fin, header);
        int numExtra = countFields(header) - 9;   // --ihs-detail columns (U3)
        for (int j = 0; j < fileLoci[i]; j++)
        {
            fin >> chrs[i];
            fin >> name;
            fin >> pos;
            fin >> freq;
            fin >> ihh1;
            fin >> ihh2;
            fin >> data;
            for (int k = 0; k < numExtra; k++) fin >> junk;
            fin >> normedData;
            fin >> crit;

            if (j > 0 && (chrs[i] != prevChr || pos < prevPos))
            {
                string what = (chrs[i] != prevChr) ? ("chromosome " + chrs[i] + " follows " + prevChr) : ("position " + to_string(pos) + " follows " + to_string(prevPos));
                cerr << "ERROR: " << normedfiles[i] << ": " << what << " (data line " << j + 1 << "). Window analysis (" << ARG_BPWIN << ") needs one chromosome per file, sorted by position.\n";
                flog << "ERROR: " << normedfiles[i] << ": " << what << " (data line " << j + 1 << "). Window analysis (" << ARG_BPWIN << ") needs one chromosome per file, sorted by position.\n";
                throw 1;
            }
            prevChr = chrs[i];
            prevPos = pos;
            if (pos < firstWinStart) continue;   // before --win-offset (B46)

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
        // B21: the window holding the file's last SNP is never closed by a
        // later SNP; write it here
        if (numSNPs > 0)
        {
            winStarts[i].push_back(winStart);
            nSNPs[i].push_back(numSNPs);
            fracCrit[i].push_back(double(numCrit) / double(numSNPs));
            if (numSNPs >= minSNPs && numCrit >= 0) numWindows++;
            maxAbsScore[i].push_back(maxAbs);
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
        for (int j = 0; j < (int)nSNPs[i].size(); j++)
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
    snum::sort2(allSNPsPerWindow, allFracCritPerWindow, numWindows);

    double *quantileBound = new double[numQuantiles];
    map<string, double> *topWindowBoundary = new map<string, double>[numQuantiles];
    windowQuantiles(allSNPsPerWindow, allFracCritPerWindow, numWindows, numQuantiles, quantileBound, topWindowBoundary, "");

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
        fout<<"chr\tstart\tend\tn_snps\tfrac_extreme\tperc\tscore\n";  // <start> <end> <num-snps> <fraction-of-extreme-snps> <percentile> <max-score of ihs/nsl etc.>
    
        for (int j = 0; j < (int)nSNPs[i].size(); j++)
        {
            if (nSNPs[i][j] < minSNPs || fracCrit[i][j] < 0)
            {
                fout << chrs[i] << "\t" << winStarts[i][j] << "\t" << winStarts[i][j] + winSize - 1 << "\t" << nSNPs[i][j] << "\t" << fracCrit[i][j] << "\t-1\tNA" << endl;
                continue;
            }
            double percentile = windowPercentile(nSNPs[i][j], fracCrit[i][j], quantileBound, numQuantiles, topWindowBoundary);

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

            fout << chrs[i] << "\t" << winStarts[i][j] << "\t" << winStarts[i][j] + winSize - 1 << "\t" << nSNPs[i][j] << "\t" << fracCrit[i][j] << "\t" << percentile << "\t";
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
    vector<string> chrs(nfiles);
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

    char str[16];
    snprintf(str, sizeof str, "%d", winSize / 1000);

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
        int winStart = 1 + this->WINOFFSET;
        int winEnd = winStart + winSize - 1;
        // B36, B46: the window analysis needs one chromosome per file in
        // position order, and starts at --win-offset
        const int firstWinStart = winStart;
        string prevChr;
        int prevPos = 0;
        int numSNPs = 0;
        int numCritTop = 0;
        int numCritBot = 0;
        double max = -99999;
        double min = 99999;
        for (int j = 0; j < fileLoci[i]; j++)
        {
            fin >> chrs[i];
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

            if (j > 0 && (chrs[i] != prevChr || pos < prevPos))
            {
                string what = (chrs[i] != prevChr) ? ("chromosome " + chrs[i] + " follows " + prevChr) : ("position " + to_string(pos) + " follows " + to_string(prevPos));
                cerr << "ERROR: " << normedfiles[i] << ": " << what << " (data line " << j + 1 << "). Window analysis (" << ARG_BPWIN << ") needs one chromosome per file, sorted by position.\n";
                flog << "ERROR: " << normedfiles[i] << ": " << what << " (data line " << j + 1 << "). Window analysis (" << ARG_BPWIN << ") needs one chromosome per file, sorted by position.\n";
                throw 1;
            }
            prevChr = chrs[i];
            prevPos = pos;
            if (pos < firstWinStart) continue;   // before --win-offset (B46)

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
        // B21: the window holding the file's last SNP is never closed by a
        // later SNP; write it here
        if (numSNPs > 0)
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
        for (int j = 0; j < (int)nSNPs[i].size(); j++)
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
    snum::sort2(allSNPsPerWindowTop, allFracCritPerWindowTop, numWindowsTop);
    snum::sort2(allSNPsPerWindowBot, allFracCritPerWindowBot, numWindowsBot);

    double *quantileBoundTop = new double[numQuantiles];
    double *quantileBoundBot = new double[numQuantiles];
    map<string, double> *topWindowBoundaryTop = new map<string, double>[numQuantiles];
    map<string, double> *topWindowBoundaryBot = new map<string, double>[numQuantiles];
    windowQuantiles(allSNPsPerWindowTop, allFracCritPerWindowTop, numWindowsTop, numQuantiles, quantileBoundTop, topWindowBoundaryTop, " (high scores)");
    windowQuantiles(allSNPsPerWindowBot, allFracCritPerWindowBot, numWindowsBot, numQuantiles, quantileBoundBot, topWindowBoundaryBot, " (low scores)");

    delete [] allSNPsPerWindowTop;
    delete [] allFracCritPerWindowTop;
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
        fout<<"chr\tstart\tend\tn_snps\tfrac_top\tfrac_bottom\tperc_top\tperc_bottom\ttop_score\tbottom_score\n";  
        cerr << "Creating window file " << winfilename[i] << endl; 
        flog << "Creating window file " << winfilename[i] << endl;

        for (int j = 0; j < (int)nSNPs[i].size(); j++)
        {
            fout << chrs[i] << "\t" << winStarts[i][j] << "\t" << winStarts[i][j] + winSize - 1 << "\t" << nSNPs[i][j] << "\t" << fracCritTop[i][j] << "\t" << fracCritBot[i][j] << "\t";
            if (nSNPs[i][j] < minSNPs)
            {
                fout << "-1\t-1\tNA\tNA" << endl;
                continue;
            }

            double percentile = windowPercentile(nSNPs[i][j], fracCritTop[i][j], quantileBoundTop, numQuantiles, topWindowBoundaryTop);
            fout << percentile << "\t";
            percentile = windowPercentile(nSNPs[i][j], fracCritBot[i][j], quantileBoundBot, numQuantiles, topWindowBoundaryBot);
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
    vector<string> chrs(nfiles);
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

    char str[16];
    snprintf(str, sizeof str, "%d", winSize / 1000);

    string name, header;
    int pos;
    double freq1, data, normedData;
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
        int winStart = 1 + this->WINOFFSET;
        int winEnd = winStart + winSize - 1;
        // B36, B46: the window analysis needs one chromosome per file in
        // position order, and starts at --win-offset
        const int firstWinStart = winStart;
        string prevChr;
        int prevPos = 0;
        int numSNPs = 0;
        int numCrit = 0;
        double maxAbs = -99999.9;   // per file (B45: carried over from the previous file)
        for (int j = 0; j < fileLoci[i]; j++)
        {
            fin >> chrs[i];
            fin >> name;
            fin >> pos;
            fin >> freq1;
            fin >> data;
            fin >> normedData;
            fin >> crit;

            if (j > 0 && (chrs[i] != prevChr || pos < prevPos))
            {
                string what = (chrs[i] != prevChr) ? ("chromosome " + chrs[i] + " follows " + prevChr) : ("position " + to_string(pos) + " follows " + to_string(prevPos));
                cerr << "ERROR: " << normedfiles[i] << ": " << what << " (data line " << j + 1 << "). Window analysis (" << ARG_BPWIN << ") needs one chromosome per file, sorted by position.\n";
                flog << "ERROR: " << normedfiles[i] << ": " << what << " (data line " << j + 1 << "). Window analysis (" << ARG_BPWIN << ") needs one chromosome per file, sorted by position.\n";
                throw 1;
            }
            prevChr = chrs[i];
            prevPos = pos;
            if (pos < firstWinStart) continue;   // before --win-offset (B46)

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
        // B21: the window holding the file's last SNP is never closed by a
        // later SNP; write it here
        if (numSNPs > 0)
        {
            winStarts[i].push_back(winStart);
            nSNPs[i].push_back(numSNPs);
            fracCrit[i].push_back(double(numCrit) / double(numSNPs));
            if (numSNPs >= minSNPs && numCrit >= 0) numWindows++;
            maxAbsScore[i].push_back(maxAbs);
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
        for (int j = 0; j < (int)nSNPs[i].size(); j++)
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
    snum::sort2(allSNPsPerWindow, allFracCritPerWindow, numWindows);

    double *quantileBound = new double[numQuantiles];
    map<string, double> *topWindowBoundary = new map<string, double>[numQuantiles];
    windowQuantiles(allSNPsPerWindow, allFracCritPerWindow, numWindows, numQuantiles, quantileBound, topWindowBoundary, "");

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
        fout<<"chr\tstart\tend\tn_snps\tfrac_extreme\tperc\tscore\n";  // <start> <end> <num-snps> <fraction-of-extreme-snps> <percentile> <max-score of ihs/nsl etc.>
    
        cerr << "Creating window file " << winfilename[i] << endl;
        flog << "Creating window file " << winfilename[i] << endl;
        for (int j = 0; j < (int)nSNPs[i].size(); j++)
        {
            if (nSNPs[i][j] < minSNPs || fracCrit[i][j] < 0)
            {
                fout << chrs[i] << "\t" << winStarts[i][j] << "\t" << winStarts[i][j] + winSize - 1 << "\t" << nSNPs[i][j] << "\t" << fracCrit[i][j] << "\t-1" <<"\tNA"<< endl;
                continue;
            }
            double percentile = windowPercentile(nSNPs[i][j], fracCrit[i][j], quantileBound, numQuantiles, topWindowBoundary);

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
            fout << chrs[i] << "\t" << winStarts[i][j] << "\t" << winStarts[i][j] + winSize - 1 << "\t" << nSNPs[i][j] << "\t" << fracCrit[i][j] << "\t" << percentile << "\t";
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
// Frequency bin of a locus, or -1 if it has none.
// Bins are [0, t0), [t0, t1), ... as before. A frequency of 1.0 is at or past
// the last edge and, before this helper, matched no bin at all: the locus was
// left out of the bin statistics and then written with the PREVIOUS row's
// normalised score and crit flag (normedData/numInBin are not reset per row).
// XP-EHH/XP-nSL/iHH12 normalise in one genome-wide bin, where p1 = 1.0 is a
// sweep completed in population 1, so for them the last bin is closed at 1.0
// (includeFreqOne). iHS/nSL bin by derived frequency, and a site fixed for the
// derived allele has no frequency class; those loci return -1, are excluded,
// and are counted.
static inline int freqBin(double f, int numBins, const double threshold[], bool includeFreqOne)
{
    if (f >= 1.0) return includeFreqOne ? numBins - 1 : -1;
    // B34: bins are left-closed [b/numBins, (b+1)/numBins). The index is
    // computed directly rather than by comparing with (b+1)*step, which in
    // floating point lands above some decimal edges (7 of 19 at 20 bins) and
    // sent frequencies printed exactly on those edges to the bin below. The
    // tolerance of 1e-6 in frequency covers selscan's 6-significant-digit
    // printing of frequencies.
    (void)threshold;
    int b = (int)std::floor(f * numBins + 1e-6 * numBins);
    if (b < 0) b = 0;
    if (b >= numBins) b = numBins - 1;
    return b;
}

// --crit-percent cutoffs (B19, U10). The quantiles are taken over every
// locus with a finite normalised score: -9999 (missing), excluded freq = 1
// loci, NaN from bins with fewer than 2 loci and loci left un-normalised by
// --log-input are left out (3.0.1 sorted all of them together, and NaN makes
// the sort order undefined). This includes loci in bins with fewer than 20
// loci, which are not written. iHH12 is one-tailed: the top critPercent is
// flagged. The other statistics flag critPercent/2 in each tail.
void SelscanNorm::critCutoffs(const double score[], int nloci, double critPercent, double critNum,
                              bool oneTailed, double &upperCutoff, double &lowerCutoff)
{
    if (!(critPercent != DEFAULT_CRIT_PERCENT && (critPercent > 0 && critPercent < 1))) {
        upperCutoff = critNum;
        lowerCutoff = -critNum;
        return;
    }
    vector<double> v;
    v.reserve(nloci);
    for (int i = 0; i < nloci; i++)
        if (std::isfinite(score[i]) && score[i] != MISSING) v.push_back(score[i]);
    if (v.empty()) {
        cerr << "ERROR: no normalised scores to compute --crit-percent cutoffs from.\n";
        flog << "ERROR: no normalised scores to compute --crit-percent cutoffs from.\n";
        exit(EXIT_FAILURE);
    }
    std::sort(v.begin(), v.end());
    if (oneTailed) {
        upperCutoff = snum::quantile_sorted(v.data(), v.size(), 1 - critPercent);
        lowerCutoff = -std::numeric_limits<double>::infinity();
    } else {
        upperCutoff = snum::quantile_sorted(v.data(), v.size(), 1 - critPercent / 2.0);
        lowerCutoff = snum::quantile_sorted(v.data(), v.size(), critPercent / 2.0);
    }
    cerr << "\n--crit-percent cutoffs from " << v.size() << " of " << nloci << " loci with a normalised score"
         << (oneTailed ? " (one-tailed)" : "") << ".\n";
    flog << "\n--crit-percent cutoffs from " << v.size() << " of " << nloci << " loci with a normalised score"
         << (oneTailed ? " (one-tailed)" : "") << ".\n";
    cerr << "Top cutoff: " << upperCutoff << endl;
    flog << "Top cutoff: " << upperCutoff << endl;
    if (!oneTailed) {
        cerr << "Bottom cutoff: " << lowerCutoff << "\n";
        flog << "Bottom cutoff: " << lowerCutoff << "\n";
    }
    cerr << "\n";
    flog << "\n";
}

int SelscanNorm::getMeanVarBins(double freq[], double data[], int nloci, double mean[], double variance[], int n[], int numBins, double threshold[], bool includeFreqOne)
{
    int excluded = 0;
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
        int b = freqBin(freq[i], numBins, threshold, includeFreqOne);
        if (b < 0) { excluded++; data[i] = MISSING; continue; }
        n[b]++;
        mean[b] += data[i];
        variance[b] += data[i] * data[i];
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
        int b = freqBin(freq[i], numBins, threshold, includeFreqOne);
        data[i] = (data[i] - mean[b]) / sqrt(variance[b]);
    }
    return excluded;
}

// Number of rows in the frequency-bin table of a norm log: the lines after
// the "bin\tnum\tmean\tvariance" header up to the first line that is not
// four fields starting with a number (B20).
int SelscanNorm::countLogBins(const std::string &binFile)
{
    std::ifstream in(binFile);
    if (!in) {
        std::cerr << "Error: cannot open bin data file " << binFile << "\n";
        exit(1);
    }
    std::string line;
    bool found = false;
    while (std::getline(in, line))
        if (line.rfind("bin\tnum\tmean\tvariance", 0) == 0) { found = true; break; }
    if (!found) return 0;
    int nb = 0;
    while (std::getline(in, line)) {
        std::istringstream iss(line);
        string a, b, c, d;
        if (!(iss >> a >> b >> c >> d) || std::isnan(parseNum(a))) break;
        nb++;
    }
    return nb;
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
                if (line.rfind("bin\tnum\tmean\tvariance", 0) == 0) {  // starts with "bin"
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
        // The table is the numBins rows right after the header, in order;
        // row b must carry the upper edge (b+1)/numBins. It ends at the
        // first line that is not "edge n mean variance" (later sections of
        // the log, such as window tables, are not read). (B20)
        int b = 0;
        while (std::getline(in, line)) {
            std::istringstream iss(line);
            string tb, tn, tm, tv;
            if (!(iss >> tb >> tn >> tm >> tv)) break;
            double bin = parseNum(tb);
            if (std::isnan(bin)) break;
            if (b >= numBins) break;
            if (std::abs(bin - threshold[b]) > 1e-6) {
                std::cerr << "ERROR: bin " << b + 1 << " of " << binFile << " has upper edge " << tb
                          << ", expected " << threshold[b] << ".\n";
                exit(1);
            }
            n[b] = atoi(tn.c_str());
            mean[b] = (tm == "nan" || tm == "-nan") ? NAN : atof(tm.c_str());
            variance[b] = (tv == "nan" || tv == "-nan") ? NAN : atof(tv.c_str());
            b++;
        }
        if (b != numBins) {
            std::cerr << "ERROR: " << binFile << " lists " << b << " frequency bins; expected " << numBins << ".\n";
            exit(1);
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
        int b = freqBin(freq[i], numBins, threshold, XPORSOFT);
        if (b < 0) { data[i] = MISSING; continue; }
        if (n[b] > 1 && !std::isnan(mean[b]) && !std::isnan(variance[b])) {
            data[i] = (data[i] - mean[b]) / std::sqrt(variance[b]);
        } else {
            data[i] = NAN;   // no statistics for this bin: not normalised (B19)
        }
    }
}

//Reads a file, calculates the normalized score, and
//outputs the original row plus normed score
void SelscanNorm::normalizeIHSDataByBins(string &filename, string &outfilename, int &fileLoci, double mean[], double variance[], int n[], int numBins, double threshold[], double upperCutoff, double lowerCutoff, bool NSL)
{
    NormTokenReader fin(filename);   // P6
    NormOut fout;
    if (!fout.open(outfilename))
    {
        cerr << "ERROR: " << outfilename << " " << strerror(errno);
        throw 1;
    }

    string name, header;
    long pos;
    string chr;
    double freq, data, normedData, ihh1, ihh2;
    int numInBin = 0;
    fin.getline(header);

    // header checked in checkStatFile before any output is opened (U3)
    if (!header.empty() && header.back() == '\r') header.pop_back();
    fout.put(header + (NSL ? "\tnorm_nsl\tcrit\n" : "\tnorm_ihs\tcrit\n"));

    // --ihs-detail adds 4 columns after the score; they are copied to the
    // .norm rows so that rows match the header (U3)
    int numExtra = countFields(header) - 7;
    std::vector<std::string> extra(numExtra);

    for (int j = 0; j < fileLoci; j++)
    {
        fin.str(chr);
        fin.str(name);
        pos = fin.integer();
        freq = fin.num();
        ihh1 = fin.num();
        ihh2 = fin.num();
        data = fin.num();   // non-finite counted in readAllIHS (B18)
        if (std::isnan(data)) data = MISSING;

        for (auto &x : extra) fin.str(x);

        if (data == MISSING || std::isnan(freq)) continue;
        {
            int b = freqBin(freq, numBins, threshold, false);
            if (b < 0) continue; // freq 1.0: excluded and counted in getMeanVarBins
            normedData = (data - mean[b]) / sqrt(variance[b]);
            numInBin = n[b];
        }

        if (numInBin >= 20)
        {
            fout.put(chr); fout.put('\t');
            fout.put(name); fout.put('\t');
            fout.putInt(pos); fout.put('\t');
            fout.putG(freq); fout.put('\t');
            fout.putG(ihh1); fout.put('\t');
            fout.putG(ihh2); fout.put('\t');
            fout.putG(data); fout.put('\t');
            for (auto &x : extra) { fout.put(x); fout.put('\t'); }
            fout.putG(normedData); fout.put('\t');
            if (normedData >= upperCutoff || normedData <= lowerCutoff) fout.put("1\n", 2);
            else fout.put("0\n", 2);
        }
    }

    if (!fout.close())
    {
        cerr << "ERROR: writing " << outfilename << " failed: " << strerror(errno) << "\n";
        throw 1;
    }
}

void SelscanNorm::normalizeXPEHHDataByBins(string &filename, string &outfilename, int &fileLoci, double mean[], double variance[], int n[], int numBins, double threshold[], double upperCutoff, double lowerCutoff, bool XPNSL)
{
    NormTokenReader fin(filename);   // P6
    NormOut fout;
    if (!fout.open(outfilename))
    {
        cerr << "ERROR: " << outfilename << " " << strerror(errno);
        throw 1;
    }

    string name, header;
    long pos;
    string chr;
    double gpos, freq1, freq2, data, normedData, ihh1, ihh2;
    int numInBin = 0;

    fin.getline(header);
    if (!header.empty() && header.back() == '\r') header.pop_back();
    if (XPNSL){
        fout.put(header + "\tnorm_xpnsl\tcrit\n");
    } else{
        fout.put(header + "\tnorm_xpehh\tcrit\n");
    }

    for (int j = 0; j < fileLoci; j++)
    {
        fin.str(chr);
        fin.str(name);
        pos = fin.integer();
        gpos = fin.num();
        freq1 = fin.num();
        ihh1 = fin.num();
        freq2 = fin.num();
        ihh2 = fin.num();
        data = fin.num();   // non-finite counted in readAllXPEHH (B18)
        if (std::isnan(data)) data = MISSING;

        if (data == MISSING || std::isnan(freq1)) continue;
        {
            int b = freqBin(freq1, numBins, threshold, true); // single bin, closed at 1.0
            normedData = (data - mean[b]) / sqrt(variance[b]);
            numInBin = n[b];
        }

        if (numInBin >= 20)
        {
            fout.put(chr); fout.put('\t');
            fout.put(name); fout.put('\t');
            fout.putInt(pos); fout.put('\t');
            fout.putG(gpos); fout.put('\t');
            fout.putG(freq1); fout.put('\t');
            fout.putG(ihh1); fout.put('\t');
            fout.putG(freq2); fout.put('\t');
            fout.putG(ihh2); fout.put('\t');
            fout.putG(data); fout.put('\t');
            fout.putG(normedData); fout.put('\t');
            if (normedData >= upperCutoff) fout.put("1\n", 2);
            else if (normedData <= lowerCutoff) fout.put("-1\n", 3);
            else fout.put("0\n", 2);
        }
    }

    if (!fout.close())
    {
        cerr << "ERROR: writing " << outfilename << " failed: " << strerror(errno) << "\n";
        throw 1;
    }
}


void SelscanNorm::normalizeIHH12DataByBins(string &filename, string &outfilename, int &fileLoci, double mean[], double variance[], int n[], int numBins, double threshold[], double upperCutoff, double /*lowerCutoff: iHH12 is one-tailed (U10)*/)
{
    NormTokenReader fin(filename);   // P6
    NormOut fout;
    if (!fout.open(outfilename))
    {
        cerr << "ERROR: " << outfilename << " " << strerror(errno);
        throw 1;
    }

    string chr, name, header;
    long pos;
    double freq1, data, normedData;
    int numInBin = 0;

    fin.getline(header);

    // header checked in checkStatFile before any output is opened (U3)
    if (!header.empty() && header.back() == '\r') header.pop_back();
    fout.put(header + "\tnorm_ihh12\tcrit\n");

    for (int j = 0; j < fileLoci; j++)
    {
        fin.str(chr);
        fin.str(name);
        pos = fin.integer();
        freq1 = fin.num();
        data = fin.num();   // non-finite counted in readAllIHH12 (B18)
        if (std::isnan(data)) data = MISSING;

        if (data == MISSING || std::isnan(freq1)) continue;
        {
            int b = freqBin(freq1, numBins, threshold, true); // single bin, closed at 1.0
            normedData = (data - mean[b]) / sqrt(variance[b]);
            numInBin = n[b];
        }

        if (numInBin >= 20)
        {
            fout.put(chr); fout.put('\t');
            fout.put(name); fout.put('\t');
            fout.putInt(pos); fout.put('\t');
            fout.putG(freq1); fout.put('\t');
            fout.putG(data); fout.put('\t');
            fout.putG(normedData); fout.put('\t');
            if (normedData >= upperCutoff) fout.put("1\n", 2);   // one-tailed (U10)
            else fout.put("0\n", 2);
        }
    }

    if (!fout.close())
    {
        cerr << "ERROR: writing " << outfilename << " failed: " << strerror(errno) << "\n";
        throw 1;
    }
}

// The exact header selscan writes for each statistic (with norm and crit
// appended for .norm files). --ihs-detail adds four iHH columns. (U3)
std::vector<std::string> SelscanNorm::expectedHeaders(const std::string &stat, bool norm)
{
    std::vector<std::string> h;
    if (stat == "ihs" || stat == "nsl") {
        std::string base = "chr\tid\tpos\tfreq\tihh1\tihh0\t" + stat;
        h.push_back(base);
        h.push_back(base + "\tihh1_left\tihh1_right\tihh0_left\tihh0_right");
    }
    else if (stat == "xpehh") h.push_back("chr\tid\tpos\tgpos\tp1\tihh1\tp2\tihh2\txpehh");
    else if (stat == "xpnsl") h.push_back("chr\tid\tpos\tgpos\tp1\tsL1\tp2\tsL2\txpnsl");
    else h.push_back("chr\tid\tpos\tfreq\tihh12");
    if (norm) for (auto &s : h) s += "\tnorm_" + stat + "\tcrit";
    return h;
}

// Checks the header and that every row has as many columns as the header,
// before any output is written (U3). Returns the number of data rows.
int SelscanNorm::checkStatFile(std::ifstream &fin, const std::string &filename, const std::string &stat, bool norm)
{
    std::streampos start = fin.tellg();
    std::string header, line;
    if (!std::getline(fin, header)) {
        cerr << "ERROR: " << filename << " is empty.\n";
        flog << "ERROR: " << filename << " is empty.\n";
        throw 0;
    }
    if (!header.empty() && header.back() == '\r') header.pop_back();

    std::vector<std::string> ok = expectedHeaders(stat, norm);
    if (std::find(ok.begin(), ok.end(), header) == ok.end()) {
        string kind = string(norm ? "normalized " : "") + "--" + stat + " files";
        cerr << "ERROR: " << filename << " has the header '" << header << "', but " << kind << " have '" << ok[0] << "'";
        flog << "ERROR: " << filename << " has the header '" << header << "', but " << kind << " have '" << ok[0] << "'";
        if (ok.size() > 1) { cerr << " or '" << ok[1] << "'"; flog << " or '" << ok[1] << "'"; }
        cerr << ".\n";
        flog << ".\n";
        if (header.compare(0, 4, "chr\t") != 0) {   // M11
            cerr << "Output from selscan before v3 has no header and no chr column. Rerun selscan 3 to make files that norm accepts.\n";
            flog << "Output from selscan before v3 has no header and no chr column. Rerun selscan 3 to make files that norm accepts.\n";
        }
        throw 0;
    }

    int ncols = countFields(header);
    int nloci = 0, lineNum = 1;
    while (std::getline(fin, line))
    {
        lineNum++;
        if (line.empty() || line == "\r") continue;
        nloci++;
        int current_cols = countFields(line);
        if (current_cols != ncols)
        {
            cerr << "ERROR: " << filename << " line " << lineNum << " has " << current_cols << " columns, but the header has " << ncols << ".\n";
            flog << "ERROR: " << filename << " line " << lineNum << " has " << current_cols << " columns, but the header has " << ncols << ".\n";
            throw 0;
        }
    }
    if (nloci == 0) {   // U11
        cerr << "WARNING: " << filename << " has a header but no data rows.\n";
        flog << "WARNING: " << filename << " has a header but no data rows.\n";
    }

    fin.clear();
    fin.seekg(start);
    return nloci;
}

void SelscanNorm::readAllIHS(vector<string> filename, int fileLoci[], int nfiles, double freq[], double score[])
{
    //chr	id	pos	freq	ihh1	ihh0	ihs
    string junk;
    int overallCount = 0;
    for (int i = 0; i < nfiles; i++)
    {
        NormTokenReader fin(filename[i]);   // P6
        fin.getline(junk); // header
        int numColsToSkip = countFields(junk) - 7; // --ihs-detail columns

        for (int j = 0; j < fileLoci[i]; j++)
        {
            fin.skip(3);
            freq[overallCount] = fin.num();
            fin.skip(2);
            double s = fin.num();
            if (std::isnan(s)) { s = MISSING; nBadRows++; }   // B18
            score[overallCount] = s;
            if (std::isnan(freq[overallCount]) && score[overallCount] != MISSING) { score[overallCount] = MISSING; nBadRows++; }
            fin.skip(numColsToSkip);
            overallCount++;
        }
    }
}

void SelscanNorm::readAllXPEHH(vector<string> filename, int fileLoci[], int nfiles, double freq1[], double freq2[], double score[])
{
    string junk;
    int overallCount = 0;
    for (int i = 0; i < nfiles; i++)
    {
        NormTokenReader fin(filename[i]);   // P6
        fin.getline(junk); // header
        for (int j = 0; j < fileLoci[i]; j++)
        {
            fin.skip(4);                    // chr id pos gpos
            freq1[overallCount] = fin.num();
            fin.skip(1);                    // ihh1
            freq2[overallCount] = fin.num();
            fin.skip(1);                    // ihh2
            double s = fin.num();
            if (std::isnan(s)) { s = MISSING; nBadRows++; }   // B18
            score[overallCount] = s;
            if (std::isnan(freq1[overallCount]) && score[overallCount] != MISSING) { score[overallCount] = MISSING; nBadRows++; }
            overallCount++;
        }
    }
}

void SelscanNorm::readAllIHH12(vector<string> filename, int fileLoci[], int nfiles, double freq1[], double score[])
{
    string junk;
    int overallCount = 0;
    for (int i = 0; i < nfiles; i++)
    {
        NormTokenReader fin(filename[i]);   // P6
        fin.getline(junk); // header
        for (int j = 0; j < fileLoci[i]; j++)
        {
            fin.skip(3);                    // chr id pos
            freq1[overallCount] = fin.num(); // p1
            double s = fin.num();            // ihh12
            if (std::isnan(s)) { s = MISSING; nBadRows++; }   // B18
            score[overallCount] = s;
            if (std::isnan(freq1[overallCount]) && score[overallCount] != MISSING) { score[overallCount] = MISSING; nBadRows++; }
            overallCount++;
        }
    }
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

