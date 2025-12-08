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


//


#include <iostream>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <algorithm> // for std::shuffle
#include <numeric>   // for std::accumulate
#include <ctime>     // for seeding



using namespace std;




struct Gene {
    std::string chrom;
    int start;
    int end;
    std::string name;
};

struct GeneNoChr {
    int start;
    int end;
    std::string name;
};



struct GeneInfo {
    std::string chrom;
    int merged_length;
    double max_score;
};

struct GeneTableEntry {
    double maxScore;
    int exonicLength = 0;
    int nWin = 0;
    int geneEnd;
};

struct Window {
    string chrom;
    int start;
    int end;
    int nSNPs;
    double frac_max;
    double frac_min;
    double perc_top;
    double perc_bottom;
    double score_min;
    double score_max;
    std::string overlap_genes;
};
// Merge intervals for a vector of intervals and return total length
int mergeIntervalsLength(std::vector<std::pair<int,int>>& intervals){
    if(intervals.empty()) return 0;
    std::sort(intervals.begin(), intervals.end());
    int total = 0;
    int curStart = intervals[0].first;
    int curEnd   = intervals[0].second;

    for(size_t i = 1; i < intervals.size(); i++){
        if(intervals[i].first > curEnd){
            total += curEnd - curStart + 1;
            curStart = intervals[i].first;
            curEnd   = intervals[i].second;
        } else {
            curEnd = std::max(curEnd, intervals[i].second);
        }
    }
    total += curEnd - curStart + 1;
    return total;
}


// added afterwards
#include <gsl/gsl_multifit.h>

struct BinData {
    std::vector<double> bins;
    std::vector<double> means;
    std::vector<double> variances;
};




class SelscanNorm{
    public:
    const string VERSION = "3.0";

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
    const string DEFAULT_FILES = "__infile__";
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

    // const string ARG_SNPWIN = "--snp-win";
    // const bool DEFAULT_SNPWIN = false;
    // const string HELP_SNPWIN = "<not implemented> If set, will use windows of a constant\n\
    // \tSNP size with varying bp length.";

    // const string ARG_SNPWINSIZE = "--snp-win-size";
    // const int DEFAULT_SNPWINSIZE = 50;
    // const string HELP_SNPWINSIZE = "<not implemented> The number of SNPs in a window.";

    const string ARG_BPWIN = "--bp-win";
    const bool DEFAULT_BPWIN = false;
    //const string HELP_BPWIN = "If set, will use windows of a constant bp size with varying\n\
    //\tnumber of SNPs.";
    const string HELP_BPWIN = "Use fixed-size bp windows (variable SNP count); outputs .windows with max/min score.\n";

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
    const string DEFAULT_LOG_INPUT = "__logfile__";
    const string HELP_LOG_INPUT =     "Specifies the log file used as an input for normalization.\n"
    "If provided, frequency-bin or mean/variance normalization is applied from the log-input file.\n"
    "Cannot be used together with --bins.\n"
    "Default: __logfile__.\n";
    const string ARG_BED = "--gene-bed"; // CHR START END GENE 
    const string DEFAULT_BED = "__filebed__"; 
    const string HELP_BED = "Provide a .bed file (<chr> <start> <end> <gene>) with gene annotations. "
                             "Used in conjunction with --annotate-win.";

    const string ARG_WIN_FILE = "--annotate-win";
    const string DEFAULT_WIN_FILE = "__filewin__";
    const string HELP_WIN_FILE = "Provide a .windows file generated by selscan norm to be annotated "
                                      "with gene names (requires --gene-bed).";

    const string ARG_GENE_SETA = "--gene-target";
    const string DEFAULT_GENE_SETA = "__genetarget__";
    const string HELP_GENE_SETA = "Provide a .genetable for target gene set for permutation tests.";

    const string ARG_GENE_SETB = "--gene-background";
    const string DEFAULT_GENE_SETB = "__genebackground__";
    const string HELP_GENE_SETB = "Provide a .genetable for background gene set for permutation test.";

    const string ARG_FINE_PERCENTILE = "--fine-percentile";
    const bool DEFAULT_FINE_PERCENTILE = false;
    const string HELP_FINE_PERCENTILE = "If set, will use fine grain percentiles (1,2,3,...,100) for normalization.";

    // const string ARG_NO_HEADER = "--no-header";
    // const bool DEFAULT_NO_HEADER = false;
    // const string HELP_NO_HEADER = "If set, will not include header in output files.";

    const int MISSING = -9999;
    const int MISSING_SCORE = 99999;

    bool FINE_PERCENTILE = false;
    string GENE_BED = "";

    bool USE_GENE_BED = false;

    //returns number of lines in file
    //throws 0 if the file fails
    int checkIHSfile(ifstream &fin);
    int checkXPEHHfile(ifstream &fin);
    int checkIHH12file(ifstream &fin);

    void readAllIHS(vector<string> filename, int fileLoci[], int nfiles, double freq[], double score[], string& chr);
    void readAllXPEHH(vector<string> filename, int fileLoci[], int nfiles, double freq1[], double freq2[], double score[], string& chr);
    void readAllIHH12(vector<string> filename, int fileLoci[], int nfiles, double freq1[], double score[], string& chr);

    void getMeanVarBins(double freq[], double data[], int nloci, double mean[], double variance[], int n[], int numBins, double threshold[]);

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
    std::vector<double> regress_out_length(const std::vector<double> &lengths,
                                       const std::vector<double> &scores);//Score_length_corrected


    std::pair<double,double> fit_length_regression(
    const std::vector<double> &lengths,
    const std::vector<double> &scores);

    std::vector<double> compute_length_residuals(
    const std::vector<double> &lengths,
    const std::vector<double> &scores,
    double beta0,
    double beta1);

    void getMeanVarBinsFromLog(const std::string &binFile,
                                       double freq[], double data[], int nloci,
                                       double mean[], double variance[], int n[],
                                       int numBins, double threshold[], bool XPORSOFT);
    
        
    double mean(const std::vector<double> &v) {
        return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
    }

        // Find column index by name
    int find_col_index(const std::vector<std::string>& headers, const std::string& col_name) {
        for (size_t i = 0; i < headers.size(); ++i)
            if (headers[i] == col_name) return i;
        return -1;
    }

    // Read file with header: returns scores vector
    std::vector<double> read_scores(const std::string& filename, const std::string& gene_col_name, const std::string& score_col_name) {
        std::ifstream f(filename);
        if (!f) {
            std::cerr << "Error opening file: " << filename << "\n";
            exit(1);
        }

        std::string line;
        std::getline(f, line); // header
        std::istringstream hss(line);
        std::vector<std::string> headers;
        std::string col;
        while (hss >> col) headers.push_back(col);

        int gene_col = find_col_index(headers, gene_col_name);
        int score_col = find_col_index(headers, score_col_name);
        if (gene_col == -1 || score_col == -1) {
            std::cerr << "Error: Columns not found in " << filename << "\n";
            exit(1);
        }

        std::vector<double> scores;
        while (std::getline(f, line)) {
            std::istringstream iss(line);
            std::vector<std::string> fields(headers.size());
            for (size_t i = 0; i < headers.size(); ++i) iss >> fields[i];
            scores.push_back(std::stod(fields[score_col]));
        }

        return scores;
    }

    void perm_test(string fileA, string fileB) {
        std::vector<double> scoresA = read_scores(fileA, "gene", "max_lenreg");
        std::vector<double> scoresB = read_scores(fileB, "gene", "max_lenreg");

        // --- Permutation test ---
        int n_perms = 10000;
        int count = 0;
        std::vector<double> combined = scoresA;
        combined.insert(combined.end(), scoresB.begin(), scoresB.end());

        int nA = scoresA.size();
        int n_total = combined.size();
        std::vector<int> indices(n_total);
        std::iota(indices.begin(), indices.end(), 0);

        gsl_rng_env_setup();
        gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_set(r, std::time(nullptr));

        double T_obs = gsl_stats_mean(scoresA.data(), 1, scoresA.size()) -
                    gsl_stats_mean(scoresB.data(), 1, scoresB.size());

        for (int i = 0; i < n_perms; ++i) {
            gsl_ran_shuffle(r, indices.data(), n_total, sizeof(int));

            std::vector<double> permA(nA), permB(n_total - nA);
            for (int j = 0; j < nA; ++j) permA[j] = combined[indices[j]];
            for (int j = nA; j < n_total; ++j) permB[j - nA] = combined[indices[j]];

            double T_perm = gsl_stats_mean(permA.data(), 1, permA.size()) -
                            gsl_stats_mean(permB.data(), 1, permB.size());
            if (T_perm >= T_obs) count++;
        }

        double p_value = static_cast<double>(count) / n_perms;
        std::cout << "Observed mean difference (A - B): " << T_obs << "\n";
        std::cout << "Empirical p-value: " << p_value << "\n";
        std::cout << ((p_value < 0.05) ? "Set A is significantly higher than Set B.\n"
                                    : "No significant difference between Set A and Set B.\n");

        gsl_rng_free(r);
    }



    void readGenes(std::string bedfile, std::vector<Gene> &genes){ 
        std::ifstream infile(bedfile);
        if (!infile) {
            std::cerr << "ERROR: could not open file " << bedfile << "\n";
            exit(EXIT_FAILURE);
        }

        std::string line;

        while (std::getline(infile, line)) {
            if (line.empty()) continue; // skip empty lines

            std::istringstream iss(line);
            std::string chrom, name;
            int start, end;

            if (!(iss >> chrom >> start >> end >> name)) {
                std::cerr << "Warning: malformed line skipped -> " << line << "\n";
                continue;
            }

            genes.push_back({chrom, start, end, name});
        }

        infile.close();

        // --- print results ---
        // for (const auto &g : genes) {
        //     std::cout << g.chrom << " " << g.name << " : "
        //             << g.start << "-" << g.end << "\n";
        // }
    }


    void annotateWindows(std::string geneBedFile, vector<std::string> windowFiles, bool XP){

        // WINDOWS integrity checks
        // DO CHECK: ALL WINDOWS in a single file MUST HAVE SAME CHR
        // DO CHECK : start < end
        // DO CHECK : consecutive line start must be > previous end

        std::ifstream geneFile(geneBedFile); //SINGLE BED FILE
        std::map<string, GeneTableEntry> geneTableMap; // AGGREGATE ALL SCORES

        // std::map<string, double> geneLength_by_gene_id; // AGGREGATE ALL LENGTHS : EXONIC GENE LENGTH : CDS LENGTH
        // //std::map<string, int> geneStart_by_gene_id; // AGGREGATE ALL LENGTHS
        // std::map<string, int> geneEnd_by_gene_id; // AGGREGATE ALL LENGTHS

//        std::map<string, std::vector<double>> geneScores; // AGGREGATE ALL SCORES: only include genes with valid scores
//        std::map<string, std::vector<double>> geneLengthsMap; // AGGREGATE ALL LENGTHS: only include genes with valid scores
        if (!geneFile) {
            std::cerr << "ERROR: could not open gene BED file " << geneBedFile << "\n";
            exit(EXIT_FAILURE);
        }
        std::vector<Gene> genes; // list of ALL genes
        
        string line;
        while (std::getline(geneFile, line)) {
            std::istringstream iss(line);
            Gene g;
            iss >> g.chrom >> g.start >> g.end >> g.name; // since BED format is 0-based start end and end exlusive
            g.end -= 1; // convert to inclusive end 
            if(g.start >= g.end){
                std::cerr << "WARNING: gene with non-positive length skipped: "<< line << "\n";
                continue;
            }
            genes.push_back(g);
        }
        // // //readGenes(GENE_BED, genes);


        std::sort(genes.begin(), genes.end(),
            [](const Gene& a, const Gene& b) {
                if (a.chrom != b.chrom) 
                    return a.chrom < b.chrom;   // lexicographic chromosome order
                if (a.start != b.start) 
                    return a.start < b.start;   // sort by start
                return a.end < b.end;           // sort by end
            }
        );

        // Organize genes by chromosome
        std::unordered_map<std::string, std::vector<GeneNoChr>> genes_by_chr;
        for(const auto& g : genes){
            GeneNoChr gnc{g.start, g.end, g.name};
            genes_by_chr[g.chrom].push_back(gnc);
        }
        //delete genes
        genes.clear();


       //std::sort(genes.begin(), genes.end(), [](auto& a, auto& b){ return a.start < b.start; }); //sorted by start position 




        for(const auto& windowFile : windowFiles){
            std::ifstream winFile(windowFile);
            if (!winFile) {
                std::cerr << "ERROR: could not open window file " << windowFile << "\n";
                exit(EXIT_FAILURE);
            }
            //open output file name appended with .annotated
            std::ofstream fout(windowFile + ".ann");
            if (!fout) {
                std::cerr << "ERROR: could not open output file " << windowFile
                            << ".ann\n";
                exit(EXIT_FAILURE);
            }
            
            std::vector<Window> windows; // all window-ranges from a single windows file
        
            std::getline(winFile, line); // skip header
            while (std::getline(winFile, line)) {
                std::istringstream iss(line);
                Window w;
                std::string score_str;
                std::string score_str_bottom;


                //fout<<"start\tend\tnSNPs\tfrac_top\tfrac_bottom\tperc\ttop_score\tbottom_score\n";  
                if(XP){                
                    //1	100001	1931	0.00517866	0	100	100	2.703	-0.623614
                    iss >> w.chrom >> w.start >> w.end >> w.nSNPs >> w.frac_max >> w.frac_min >> w.perc_top >> w.perc_bottom >> score_str >> score_str_bottom;
                }else{
                    iss >> w.chrom >> w.start >> w.end >> w.nSNPs >> w.frac_max >> w.perc_top >> score_str;
                    //cout<<"Read window: "<< w.start << "-" << w.end << " nSNPs: "<< w.nSNPs << " frac_max: "<< w.frac_max << " perc: "<< w.perc << " score: "<< score_str << "\n";
                }
                w.score_max = (score_str == "NA") ? -MISSING_SCORE : std::stod(score_str);
                if(XP){
                    w.score_min = (score_str_bottom == "NA") ? MISSING_SCORE : std::stod(score_str_bottom);
                }

                w.overlap_genes = "";
                windows.push_back(w);

                //check chromosome, only get vector of genes for that chromosome
                std::string chrom_from_win = w.chrom; // Need to get chromosome from window file if available
                auto it = genes_by_chr.find(chrom_from_win);
                if (it == genes_by_chr.end()) {
                    std::cerr << "ERROR: chromosome " << chrom_from_win << " not found in gene BED file.\n";
                    exit(EXIT_FAILURE);
                }
                const auto& genes = it->second;

                 // Two-pointer approach
                size_t gene_idx = 0;
                for (auto &w : windows) {
                    std::string ov;

                    // Advance gene_idx to the first gene that might overlap
                    while (gene_idx < genes.size() && genes[gene_idx].end <= w.start) {
                        gene_idx++;
                    }

                    // Check overlapping genes starting from current gene_idx
                    size_t j = gene_idx;
                    string GENE_ID = chrom_from_win + ":" + genes[j].name ;

                    while (j < genes.size() && genes[j].start < w.end) {
                        if ((w.end > genes[j].start && w.start < genes[j].end)) {  // overlaps 
                            if (!ov.empty()) ov += ", ";
                            ov += genes[j].name;

                            //geneScores[GENE_ID].push_back(w.score_max);
                            // geneScores[genes[j].name].push_back(w.score_max);
                            // geneLengthsMap[genes[j].name].push_back(genes[j].end - genes[j].start + 1);
                        
                            if(true){  //if((w.score_max) != MISSING_SCORE){
                                // if a region did not see any overlap (no scores assigned) we exlude that from analysis
                                if(geneTableMap.find(GENE_ID) != geneTableMap.end()){  //if exists, check if new score is higher, if higher, update
                                    double existing_score = geneTableMap[GENE_ID].maxScore;
                                    if(w.score_max > existing_score){
                                        geneTableMap[GENE_ID].maxScore = w.score_max ;
                                    }
                                    
                                    if(genes[j].start > geneTableMap[GENE_ID].geneEnd){ // non-overlapping
                                        geneTableMap[GENE_ID].exonicLength += genes[j].end - genes[j].start + 1;
                                        //geneStart_by_gene_id[GENE_ID] = genes[j].start;
                                        geneTableMap[GENE_ID].geneEnd = genes[j].end;
                                    } else { // new gene region overlaps with previous
                                        if( genes[j].end > geneTableMap[GENE_ID].geneEnd){
                                            geneTableMap[GENE_ID].exonicLength += genes[j].end - geneTableMap[GENE_ID].geneEnd;
                                            geneTableMap[GENE_ID].geneEnd = genes[j].end;
                                        }
                                    }
                                            // total += curEnd - curStart + 1;
                                }else{ //first score for that gene
                                    geneTableMap[GENE_ID].maxScore = w.score_max;
                                    geneTableMap[GENE_ID].exonicLength =  genes[j].end - genes[j].start + 1; //inlcusive both
                                    geneTableMap[GENE_ID].geneEnd = genes[j].end;
                                }
                                geneTableMap[GENE_ID].nWin += 1;
                            }

                        }
                        j++;
                    }
                    

                    w.overlap_genes = ov=="" ? "-" : ov;
                    
                }
            }

            

            //print header
            if(XP){
                fout<<"start\tend\tnSNPs\tfrac_top\tfrac_bottom\tperc\ttop_score\tbottom_score\toverlap_genes\n";  
            }else{
                fout<<"start\tend\tnSNPs\tfrac_extreme\tperc\tscore\toverlap_genes\n";  
            }
            // Print output
            for (auto &w : windows) {
                std::string score_str = (w.score_max==-1) ? "NA" : std::to_string(w.score_max); //convert score_max to string, "NA" if nan
                
                if(XP){
                    std::string score_str_bottom = (w.score_min==-1) ? "NA" : std::to_string(w.score_min); //convert score_min to string, "NA" if nan

                    fout<< w.chrom << "\t" << w.start << "\t" << w.end << "\t" << w.nSNPs << "\t"
                        << w.frac_max << "\t" << w.frac_min << "\t" << w.perc_top << "\t" << w.perc_bottom  << "\t"
                        << score_str << "\t" << score_str_bottom
                        << "\t" << w.overlap_genes << "\n";
                }else{
                    fout<< w.chrom << "\t" << w.start << "\t" << w.end << "\t" << w.nSNPs << "\t"
                        << w.frac_max << "\t" << w.perc_top << "\t"
                        << score_str
                        << "\t" << w.overlap_genes << "\n";
                }
                
            }
            fout.close();
            cout<<"Written gene annotated windows to "<< windowFile + ".ann\n";
        }


        bool PRINT_GENE_TABLE = true;
        if(PRINT_GENE_TABLE){
            vector<double> geneLengths;
            vector<double> geneScoresMax;

            for (const auto& [gene, row] : geneTableMap) { // do it in order of ID
                geneLengths.push_back(row.exonicLength); 
                geneScoresMax.push_back(row.maxScore); 
            }   
            ////ASSUMPTION BOTH ARE ORDERED         
            // for (const auto& [gene, maxscore] : geneScores_by_gene_id) { // do it in order of ID
            //     geneScoresMax.push_back(maxscore); 
            // }

            pair<double,double> beta = fit_length_regression(geneLengths, geneScoresMax);
            
            for(const auto& windowFile: windowFiles){
                cout<<"Length regression for gene max scores from file "<< windowFile <<": beta0="<< beta.first << " beta1="<< beta.second << "\n";
                ofstream genetable; // output gene table with mean, sd, nwin, max, max_lenreg
                string genetablefile = windowFile + ".genetable";
                genetable.open(genetablefile.c_str());
                if (genetable.fail())
                {
                    cerr << "ERROR: " << genetablefile << " " << strerror(errno);
                    exit(EXIT_FAILURE);
                }
                genetable << "chr\tgene\tlen\tnwin\tmax\tmax_lenreg\n";  //gene   mean   sd   nwin   max   max_lenreg
                int i = 0;

                string winchr = "chr1";
                std::string prefix = winchr+ ":";
                auto it = geneTableMap.lower_bound(prefix);  // first key >= "a:"

                while (it != geneTableMap.end() && it->first.compare(0, prefix.size(), prefix) == 0) { // while key starts with "chr1:"
                    std::cout << it->first << " " ;

                    std::string geneName = it->first.substr(prefix.size());  // remove prefix
                    GeneTableEntry entry = it->second;

                    double pred = beta.first + beta.second * std::log(entry.exonicLength);
                    double residuals = entry.maxScore - pred;

                    if(abs(entry.maxScore)!=MISSING_SCORE) genetable << winchr<< "\t" << geneName << "\t" << entry.exonicLength << "\t" << entry.nWin << "\t" << entry.maxScore << "\t" << residuals << "\t" <<"\n";
                    
                    ++it;
                }


                // for (const auto& [gene, vals] : geneScores) {
                //     double* data = const_cast<double*>(vals.data()); // gsl requires non-const pointer
                //     size_t n = vals.size();
                //     double maxScore  = *std::max_element(vals.begin(), vals.end());
                //     int len = geneLengthsMap[gene][0];
                //     // double meanScore = gsl_stats_mean(data, 1, n);
                //     // double varScore  = 0; // for n = 1
                //     // if (n > 1) varScore = std::sqrt(gsl_stats_variance(data, 1, n)/n); 
                //     //genetable << gene << "\t" << meanScore << "\t" << varScore << "\t" << n << "\t" << maxScore << "\t" << geneScoresMax[i] << "\t" <<"\n";
                //     if(abs(maxScore)!=MISSING_SCORE) genetable << gene << "\t" << len << "\t" << n << "\t" << maxScore << "\t" << geneScoresMax[i] << "\t" <<"\n";
                //     i++;
                // }
                genetable.close();
                cout<<"Written gene table to "<< windowFile + ".genetable" <<endl;
            }
        }

        //////
        // std::ifstream winFile(windowFile);
        // if (!winFile) {
        //     std::cerr << "ERROR: could not open window file " << windowFile << "\n";
        //     exit(EXIT_FAILURE);
        // }

        // //open output file name appended with .annotated
        // std::ofstream fout(windowFile + ".ann");
        // if (!fout) {
        //    std::cerr << "ERROR: could not open output file " << windowFile
        //              << ".ann\n";
        //    exit(EXIT_FAILURE);
        // }

        // std::vector<Window> windows; // all window-ranges from a single windows file
        
        // Read windows (assumes sorted by start)
        // std::getline(winFile, line); // skip header
        // while (std::getline(winFile, line)) {
        //     std::istringstream iss(line);
        //     Window w;
        //     std::string score_str;
        //     std::string score_str_bottom;


        //     //fout<<"start\tend\tnSNPs\tfrac_top\tfrac_bottom\tperc\ttop_score\tbottom_score\n";  
        //     if(XP){                
        //         //1	100001	1931	0.00517866	0	100	100	2.703	-0.623614
        //         iss >> w.start >> w.end >> w.nSNPs >> w.frac_max >> w.frac_min >> w.perc_top >> w.perc_bottom >> score_str >> score_str_bottom;
        //     }else{
        //         iss >> w.start >> w.end >> w.nSNPs >> w.frac_max >> w.perc_top >> score_str;
        //         //cout<<"Read window: "<< w.start << "-" << w.end << " nSNPs: "<< w.nSNPs << " frac_max: "<< w.frac_max << " perc: "<< w.perc << " score: "<< score_str << "\n";
        //     }
        //     w.score_max = (score_str == "NA") ? -MISSING_SCORE : std::stod(score_str);
        //     if(XP){
        //         w.score_min = (score_str_bottom == "NA") ? MISSING_SCORE : std::stod(score_str_bottom);
        //     }

        //     w.overlap_genes = "";
        //     windows.push_back(w);
        // }



        // // Two-pointer approach
        // size_t gene_idx = 0;
        // for (auto &w : windows) {
        //     std::string ov;

        //     // Advance gene_idx to the first gene that might overlap
        //     while (gene_idx < genes.size() && genes[gene_idx].end <= w.start) {
        //         gene_idx++;
        //     }

        //     // Check overlapping genes starting from current gene_idx
        //     size_t j = gene_idx;
        //     while (j < genes.size() && genes[j].start < w.end) {
        //         if ((w.end > genes[j].start && w.start < genes[j].end)) {  // overlaps 
        //             if (!ov.empty()) ov += ", ";
        //             ov += genes[j].name;

        //             // geneScores[genes[j].name].push_back(w.score_max);
        //             // geneLengthsMap[genes[j].name].push_back(genes[j].end - genes[j].start + 1);
                   
        //             if((w.score_max) != MISSING_SCORE){
        //                 if(geneScores_by_gene_id.find(genes[j].name) != geneScores_by_gene_id.end()){  //if exists, check if new score is higher, if higher, update
        //                     double existing_score = geneScores_by_gene_id[genes[j].name];
        //                     if(w.score_max > existing_score){
        //                         geneScores_by_gene_id[genes[j].name] = w.score_max ;
        //                     }
        //                 }else{
        //                     geneScores_by_gene_id[genes[j].name] = w.score_max;
        //                     geneLength_by_gene_id[genes[j].name] =  genes[j].end - genes[j].start + 1;
        //                     //update adjusted length with merged ;ength
        //                             // int total = 0;
        //                             // int curStart = genes[0].start;
        //                             // int curEnd   = genes[0].end;

        //                             // for(size_t i = 1; i < genes.size(); i++){
        //                             //     if(genes[i].start > curEnd){
        //                             //         total += curEnd - curStart + 1;
        //                             //         curStart = genes[i].start;
        //                             //         curEnd   = genes[i].end;
        //                             //     } else {
        //                             //         curEnd = std::max(curEnd, genes[i].end);
        //                             //     }
        //                             // }
        //                             // total += curEnd - curStart + 1;
        //                 }
                        
        //             }

        //         }
        //         j++;
        //     }

        //     w.overlap_genes = ov=="" ? "-" : ov;
            
        // }


        // bool PRINT_GENE_TABLE = true;
        // if(PRINT_GENE_TABLE){
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
        //     string genetablefile = windowFile + ".genetable";
        //     genetable.open(genetablefile.c_str());
        //     if (genetable.fail())
        //     {
        //         cerr << "ERROR: " << genetablefile << " " << strerror(errno);
        //         exit(EXIT_FAILURE);
        //     }
        //     genetable << "gene\tlen\tnwin\tmax\tmax_lenreg\n";  //gene   mean   sd   nwin   max   max_lenreg
        //     int i = 0;
        //     for (const auto& [gene, vals] : geneScores) {
        //         double* data = const_cast<double*>(vals.data()); // gsl requires non-const pointer
        //         size_t n = vals.size();
        //         double maxScore  = *std::max_element(vals.begin(), vals.end());
        //         int len = geneLengthsMap[gene][0];
        //         // double meanScore = gsl_stats_mean(data, 1, n);
        //         // double varScore  = 0; // for n = 1
        //         // if (n > 1) varScore = std::sqrt(gsl_stats_variance(data, 1, n)/n); 
        //         //genetable << gene << "\t" << meanScore << "\t" << varScore << "\t" << n << "\t" << maxScore << "\t" << geneScoresMax[i] << "\t" <<"\n";
        //         if(abs(maxScore)!=MISSING_SCORE) genetable << gene << "\t" << len << "\t" << n << "\t" << maxScore << "\t" << geneScoresMax[i] << "\t" <<"\n";
        //         i++;
        //     }
        //     genetable.close();
        // }

        // //print header
        // if(XP){
        //     fout<<"start\tend\tnSNPs\tfrac_top\tfrac_bottom\tperc\ttop_score\tbottom_score\toverlap_genes\n";  
        // }else{
        //     fout<<"start\tend\tnSNPs\tfrac_extreme\tperc\tscore\toverlap_genes\n";  
        // }
        // // Print output
        // for (auto &w : windows) {
        //     std::string score_str = (w.score_max==-1) ? "NA" : std::to_string(w.score_max); //convert score_max to string, "NA" if nan
            
        //     if(XP){
        //         std::string score_str_bottom = (w.score_min==-1) ? "NA" : std::to_string(w.score_min); //convert score_min to string, "NA" if nan

        //         fout<< w.start << "\t" << w.end << "\t" << w.nSNPs << "\t"
        //             << w.frac_max << "\t" << w.frac_min << "\t" << w.perc_top << "\t" << w.perc_bottom  << "\t"
        //             << score_str << "\t" << score_str_bottom
        //             << "\t" << w.overlap_genes << "\n";
        //     }else{
        //         fout<< w.start << "\t" << w.end << "\t" << w.nSNPs << "\t"
        //             << w.frac_max << "\t" << w.perc_top << "\t"
        //             << score_str
        //             << "\t" << w.overlap_genes << "\n";
        //     }
            
        // }
        // fout.close();
        // cout<<"Written gene table to "<< windowFile + ".genetable" <<endl;
        // cout<<"Written gene annotated windows to "<< windowFile + ".ann\n";
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


