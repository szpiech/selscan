#ifndef __GENE_H__
#define __GENE_H__

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

#include <iostream>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <algorithm> // for std::shuffle
#include <numeric>   // for std::accumulate
#include <ctime>     // for seeding

#include <fstream>

// added afterwards
#include <gsl/gsl_multifit.h>



using namespace std;

struct Window {
    std::string chrom;
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

struct Gene
{
    std::string chrom;
    int start;
    int end;
    std::string name;
};

struct GeneNoChr
{
    int start;
    int end;
    std::string name;
};

// struct GeneInfo
// {
//     std::string chrom;
//     int merged_length;
//     double max_score;
// };

struct GeneTableEntry
{
    double maxScore = -99999;
    int lengthSpan = 0;
    int nWin = 0;
    //int geneEnd;
    int windowFileId = -1;
};

struct GeneTableEntryTrue
{
    double maxScore = -99999;
    int lengthSpan = 0;
    int nWin = 0;
    string geneId;
};

class GeneAnalyzer
{
private:
    void readGenesFromBed(std::string bedfile, std::vector<Gene> &genes, bool canonical = false);
    void readGenesFromGTF(std::string gtffile, std::vector<Gene> &genes,  bool useTranscripts = false, bool canonical = false);
    std::vector<Gene> genes;
    const int MISSING_SCORE = 99999;   

public:

    inline double compute_length_residual(
    const double length,
    const double score,
    double beta0,
    double beta1);

    std::vector<double> regress_out_length(const std::vector<double> &lengths,
                                           const std::vector<double> &scores); // Score_length_corrected

    std::pair<double, double> fit_length_regression(
        const std::vector<double> &lengths,
        const std::vector<double> &scores);

    std::vector<double> compute_length_residuals(
        const std::vector<double> &lengths,
        const std::vector<double> &scores,
        double beta0,
        double beta1);

    void perm_test(std::string fileA, std::string fileB);

    // Merge intervals for a vector of intervals and return total length
    int mergeIntervalsLength(std::vector<std::pair<int, int>> &intervals)
    {
        if (intervals.empty())
            return 0;
        std::sort(intervals.begin(), intervals.end());
        int total = 0;
        int curStart = intervals[0].first;
        int curEnd = intervals[0].second;

        for (size_t i = 1; i < intervals.size(); i++)
        {
            if (intervals[i].first > curEnd)
            {
                total += curEnd - curStart + 1;
                curStart = intervals[i].first;
                curEnd = intervals[i].second;
            }
            else
            {
                curEnd = std::max(curEnd, intervals[i].second);
            }
        }
        total += curEnd - curStart + 1;
        return total;
    }


    // Find column index by name
    int find_col_index(const std::vector<std::string> &headers, const std::string &col_name)
    {
        for (size_t i = 0; i < headers.size(); ++i)
            if (headers[i] == col_name)
                return i;
        return -1;
    }

    // Read file with header: returns scores vector
    std::vector<double> read_scores(const std::string &filename, const std::string &gene_col_name, const std::string &score_col_name)
    {
        std::ifstream f(filename);
        if (!f)
        {
            std::cerr << "Error opening file: " << filename << "\n";
            exit(1);
        }

        std::string line;
        std::getline(f, line); // header
        std::istringstream hss(line);
        std::vector<std::string> headers;
        std::string col;
        while (hss >> col)
            headers.push_back(col);

        int gene_col = find_col_index(headers, gene_col_name);
        int score_col = find_col_index(headers, score_col_name);
        if (gene_col == -1 || score_col == -1)
        {
            std::cerr << "Error: Columns not found in " << filename << "\n";
            exit(1);
        }

        std::vector<double> scores;
        while (std::getline(f, line))
        {
            std::istringstream iss(line);
            std::vector<std::string> fields(headers.size());
            for (size_t i = 0; i < headers.size(); ++i)
                iss >> fields[i];
            scores.push_back(std::stod(fields[score_col]));
        }

        return scores;
    }

    void annotateWindows(std::string geneBedFile, vector<std::string> windowFiles, bool XP);
};

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


#endif