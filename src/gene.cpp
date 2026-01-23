
#include "gene.h"

void GeneAnalyzer::readGenes(std::string bedfile, std::vector<Gene> &genes){ 
    std::ifstream infile(bedfile);
    if (!infile) {
        std::cerr << "ERROR: could not open file " << bedfile << "\n";
        std::exit(EXIT_FAILURE);
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


void GeneAnalyzer::perm_test(string fileA, string fileB) {
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



std::vector<double> GeneAnalyzer::compute_length_residuals(
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
std::vector<double> GeneAnalyzer::regress_out_length(
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


std::pair<double,double> GeneAnalyzer::fit_length_regression(
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

    




//     // if any part of gene overlaps a window, mark that window as overlapping the gene
// // if a window overlaps multiple genes, list all genes in the annotation (comma-separated?), 
// /**
//  * Regresses score ~ intercept + log(length) using GSL
//  * and returns residuals (length-corrected scores).
//  *
//  * @param lengths  Vector of gene lengths
//  * @param scores   Vector of per-gene scores
//  * @return residuals vector (score corrected for length)
//  */
// std::vector<double> GeneAnalyzer::regress_out_length(
//     const std::vector<double> &lengths,
//     const std::vector<double> &scores) 
// {
//     size_t n = lengths.size();
//     size_t p = 2; // intercept + log(length)

//     // First collect only valid (score != -1) entries
//     std::vector<size_t> valid_idx;
//     valid_idx.reserve(n);
//     for (size_t i = 0; i < n; i++) {
//         if (abs(scores[i]) != MISSING_SCORE) {
//             valid_idx.push_back(i);
//         }
//     }

//     size_t m = valid_idx.size(); // number of valid points
//     if (m == 0) {
//         // nothing to fit, just return original
//         return scores;
//     }

//     gsl_matrix *X = gsl_matrix_alloc(m, p);
//     gsl_vector *y = gsl_vector_alloc(m);
//     gsl_vector *c = gsl_vector_alloc(p);
//     gsl_matrix *cov = gsl_matrix_alloc(p, p);
//     double chisq;

//     // Fill design matrix and response with only valid data
//     for (size_t j = 0; j < m; j++) {
//         size_t i = valid_idx[j];
//         gsl_matrix_set(X, j, 0, 1.0);                  // intercept
//         gsl_matrix_set(X, j, 1, std::log(lengths[i])); // log(length)
//         gsl_vector_set(y, j, scores[i]);
//     }

//     // Fit linear regression
//     gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(m, p);
//     gsl_multifit_linear(X, y, c, cov, &chisq, work);
//     gsl_multifit_linear_free(work);

//     double beta0 = gsl_vector_get(c, 0);
//     double beta1 = gsl_vector_get(c, 1);

//     // Compute residuals (keep same size as input)
//     std::vector<double> residuals(n, -1.0);
//     for (size_t j = 0; j < m; j++) {
//         size_t i = valid_idx[j];
//         double pred = beta0 + beta1 * std::log(lengths[i]);
//         residuals[i] = scores[i] - pred;
//     }

//     // Shift residuals so the minimum valid one is 0
//     // double minResidual = std::numeric_limits<double>::max();
//     // for (double r : residuals) {
//     //     if (r != -1.0) minResidual = std::min(minResidual, r);
//     // }
//     // if (minResidual < 0.0) {
//     //     for (auto &r : residuals) {
//     //         if (r != -1.0) r -= minResidual;
//     //     }
//     // }

//     // Free memory
//     gsl_matrix_free(X);
//     gsl_matrix_free(cov);
//     gsl_vector_free(y);
//     gsl_vector_free(c);

//     return residuals;
// }




void GeneAnalyzer::annotateWindows(std::string geneBedFile, vector<std::string> windowFiles, bool XP){

    // WINDOWS integrity checks
    // DO CHECK: ALL WINDOWS in a single file MUST HAVE SAME CHR
    // DO CHECK : start < end
    // DO CHECK : consecutive line start must be > previous end

    std::ifstream geneFile(geneBedFile);           // SINGLE BED FILE
    std::map<std::string, GeneTableEntry> geneTableMap; // AGGREGATE ALL SCORES

    // std::map<string, double> geneLength_by_gene_id; // AGGREGATE ALL LENGTHS : EXONIC GENE LENGTH : CDS LENGTH
    // //std::map<string, int> geneStart_by_gene_id; // AGGREGATE ALL LENGTHS
    // std::map<string, int> geneEnd_by_gene_id; // AGGREGATE ALL LENGTHS

    //        std::map<string, std::vector<double>> geneScores; // AGGREGATE ALL SCORES: only include genes with valid scores
    //        std::map<string, std::vector<double>> geneLengthsMap; // AGGREGATE ALL LENGTHS: only include genes with valid scores
    if (!geneFile)
    {
        std::cerr << "ERROR: could not open gene BED file " << geneBedFile << "\n";
        exit(EXIT_FAILURE);
    }
    std::vector<Gene> genes; // list of ALL genes

    std::string line;
    while (std::getline(geneFile, line))
    {
        std::istringstream iss(line);
        Gene g;
        iss >> g.chrom >> g.start >> g.end >> g.name; // since BED format is 0-based start end and end exlusive
        g.end -= 1;                                   // convert to inclusive end
        if (g.start >= g.end)
        {
            std::cerr << "WARNING: gene with non-positive length skipped: " << line << "\n";
            continue;
        }
        genes.push_back(g);
    }
    // // //readGenes(GENE_BED, genes);

    std::sort(genes.begin(), genes.end(),
                [](const Gene &a, const Gene &b)
                {
                    if (a.chrom != b.chrom)
                        return a.chrom < b.chrom; // lexicographic chromosome order
                    if (a.start != b.start)
                        return a.start < b.start; // sort by start
                    return a.end < b.end;         // sort by end
                });

    // Organize genes by chromosome
    std::unordered_map<std::string, std::vector<GeneNoChr>> genes_by_chr;
    for (const auto &g : genes)
    {
        GeneNoChr gnc{g.start, g.end, g.name};
        genes_by_chr[g.chrom].push_back(gnc);
    }
    // delete genes
    genes.clear();

    // std::sort(genes.begin(), genes.end(), [](auto& a, auto& b){ return a.start < b.start; }); //sorted by start position

    for (const auto &windowFile : windowFiles)
    {
        std::ifstream winFile(windowFile);
        if (!winFile)
        {
            std::cerr << "ERROR: could not open window file " << windowFile << "\n";
            exit(EXIT_FAILURE);
        }
        // open output file name appended with .annotated
        std::ofstream fout(windowFile + ".ann");
        if (!fout)
        {
            std::cerr << "ERROR: could not open output file " << windowFile
                        << ".ann\n";
            exit(EXIT_FAILURE);
        }

        std::vector<Window> windows; // all window-ranges from a single windows file

        std::getline(winFile, line); // skip header
        while (std::getline(winFile, line))
        {
            std::istringstream iss(line);
            Window w;
            std::string score_str;
            std::string score_str_bottom;

            // fout<<"start\tend\tnSNPs\tfrac_top\tfrac_bottom\tperc\ttop_score\tbottom_score\n";
            if (XP)
            {
                // 1	100001	1931	0.00517866	0	100	100	2.703	-0.623614
                iss >> w.chrom >> w.start >> w.end >> w.nSNPs >> w.frac_max >> w.frac_min >> w.perc_top >> w.perc_bottom >> score_str >> score_str_bottom;
            }
            else
            {
                iss >> w.chrom >> w.start >> w.end >> w.nSNPs >> w.frac_max >> w.perc_top >> score_str;
                // cout<<"Read window: "<< w.start << "-" << w.end << " nSNPs: "<< w.nSNPs << " frac_max: "<< w.frac_max << " perc: "<< w.perc << " score: "<< score_str << "\n";
            }
            w.score_max = (score_str == "NA") ? -MISSING_SCORE : std::stod(score_str);
            if (XP)
            {
                w.score_min = (score_str_bottom == "NA") ? MISSING_SCORE : std::stod(score_str_bottom);
            }

            w.overlap_genes = "";
            windows.push_back(w);

            // check chromosome, only get vector of genes for that chromosome
            std::string chrom_from_win = w.chrom; // Need to get chromosome from window file if available
            auto it = genes_by_chr.find(chrom_from_win);
            if (it == genes_by_chr.end())
            {
                std::cerr << "ERROR: chromosome " << chrom_from_win << " not found in gene BED file.\n";
                exit(EXIT_FAILURE);
            }
            const auto &genes = it->second;

            // Two-pointer approach
            size_t gene_idx = 0;
            for (auto &w : windows)
            {
                std::string ov;

                // Advance gene_idx to the first gene that might overlap
                while (gene_idx < genes.size() && genes[gene_idx].end <= w.start)
                {
                    gene_idx++;
                }

                // Check overlapping genes starting from current gene_idx
                size_t j = gene_idx;
                std::string GENE_ID = chrom_from_win + ":" + genes[j].name;

                while (j < genes.size() && genes[j].start < w.end)
                {
                    if ((w.end > genes[j].start && w.start < genes[j].end))
                    { // overlaps
                        if (!ov.empty())
                            ov += ", ";
                        ov += genes[j].name;

                        // geneScores[GENE_ID].push_back(w.score_max);
                        //  geneScores[genes[j].name].push_back(w.score_max);
                        //  geneLengthsMap[genes[j].name].push_back(genes[j].end - genes[j].start + 1);

                        if (true)
                        { // if((w.score_max) != MISSING_SCORE){
                            // if a region did not see any overlap (no scores assigned) we exlude that from analysis
                            if (geneTableMap.find(GENE_ID) != geneTableMap.end())
                            { // if exists, check if new score is higher, if higher, update
                                double existing_score = geneTableMap[GENE_ID].maxScore;
                                if (w.score_max > existing_score)
                                {
                                    geneTableMap[GENE_ID].maxScore = w.score_max;
                                }

                                if (genes[j].start > geneTableMap[GENE_ID].geneEnd)
                                { // non-overlapping
                                    geneTableMap[GENE_ID].exonicLength += genes[j].end - genes[j].start + 1;
                                    // geneStart_by_gene_id[GENE_ID] = genes[j].start;
                                    geneTableMap[GENE_ID].geneEnd = genes[j].end;
                                }
                                else
                                { // new gene region overlaps with previous
                                    if (genes[j].end > geneTableMap[GENE_ID].geneEnd)
                                    {
                                        geneTableMap[GENE_ID].exonicLength += genes[j].end - geneTableMap[GENE_ID].geneEnd;
                                        geneTableMap[GENE_ID].geneEnd = genes[j].end;
                                    }
                                }
                                // total += curEnd - curStart + 1;
                            }
                            else
                            { // first score for that gene
                                geneTableMap[GENE_ID].maxScore = w.score_max;
                                geneTableMap[GENE_ID].exonicLength = genes[j].end - genes[j].start + 1; // inlcusive both
                                geneTableMap[GENE_ID].geneEnd = genes[j].end;
                            }
                            geneTableMap[GENE_ID].nWin += 1;
                        }
                    }
                    j++;
                }

                w.overlap_genes = ov == "" ? "-" : ov;
            }
        }

        // print header
        if (XP)
        {
            fout << "start\tend\tnSNPs\tfrac_top\tfrac_bottom\tperc\ttop_score\tbottom_score\toverlap_genes\n";
        }
        else
        {
            fout << "start\tend\tnSNPs\tfrac_extreme\tperc\tscore\toverlap_genes\n";
        }
        // Print output
        for (auto &w : windows)
        {
            std::string score_str = (w.score_max == -1) ? "NA" : std::to_string(w.score_max); // convert score_max to string, "NA" if nan

            if (XP)
            {
                std::string score_str_bottom = (w.score_min == -1) ? "NA" : std::to_string(w.score_min); // convert score_min to string, "NA" if nan

                fout << w.chrom << "\t" << w.start << "\t" << w.end << "\t" << w.nSNPs << "\t"
                        << w.frac_max << "\t" << w.frac_min << "\t" << w.perc_top << "\t" << w.perc_bottom << "\t"
                        << score_str << "\t" << score_str_bottom
                        << "\t" << w.overlap_genes << "\n";
            }
            else
            {
                fout << w.chrom << "\t" << w.start << "\t" << w.end << "\t" << w.nSNPs << "\t"
                        << w.frac_max << "\t" << w.perc_top << "\t"
                        << score_str
                        << "\t" << w.overlap_genes << "\n";
            }
        }
        fout.close();
        cout << "Written gene annotated windows to " << windowFile + ".ann\n";
    }

    bool PRINT_GENE_TABLE = true;
    if (PRINT_GENE_TABLE)
    {
        vector<double> geneLengths;
        vector<double> geneScoresMax;

        for (const auto &[gene, row] : geneTableMap)
        { // do it in order of ID
            geneLengths.push_back(row.exonicLength);
            geneScoresMax.push_back(row.maxScore);
        }
        ////ASSUMPTION BOTH ARE ORDERED
        // for (const auto& [gene, maxscore] : geneScores_by_gene_id) { // do it in order of ID
        //     geneScoresMax.push_back(maxscore);
        // }

        pair<double, double> beta = fit_length_regression(geneLengths, geneScoresMax);

        for (const auto &windowFile : windowFiles)
        {
            cout << "Length regression for gene max scores from file " << windowFile << ": beta0=" << beta.first << " beta1=" << beta.second << "\n";
            ofstream genetable; // output gene table with mean, sd, nwin, max, max_lenreg
            string genetablefile = windowFile + ".genetable";
            genetable.open(genetablefile.c_str());
            if (genetable.fail())
            {
                cerr << "ERROR: " << genetablefile << " " << strerror(errno);
                exit(EXIT_FAILURE);
            }
            genetable << "chr\tgene\tlen\tnwin\tmax\tmax_lenreg\n"; // gene   mean   sd   nwin   max   max_lenreg
            int i = 0;

            string winchr = "chr1";
            std::string prefix = winchr + ":";
            auto it = geneTableMap.lower_bound(prefix); // first key >= "a:"

            while (it != geneTableMap.end() && it->first.compare(0, prefix.size(), prefix) == 0)
            { // while key starts with "chr1:"
                std::cout << it->first << " ";

                std::string geneName = it->first.substr(prefix.size()); // remove prefix
                GeneTableEntry entry = it->second;

                double pred = beta.first + beta.second * std::log(entry.exonicLength);
                double residuals = entry.maxScore - pred;

                if (abs(entry.maxScore) != MISSING_SCORE)
                    genetable << winchr << "\t" << geneName << "\t" << entry.exonicLength << "\t" << entry.nWin << "\t" << entry.maxScore << "\t" << residuals << "\t" << "\n";

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
            cout << "Written gene table to " << windowFile + ".genetable" << endl;
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