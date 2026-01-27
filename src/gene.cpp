
#include "gene.h"
#include <unordered_map>
#include "gzstream.h"

// deduplicated genes from bed, default take gene span
void GeneAnalyzer::readGenesFromBed(std::string bedfile, std::vector<Gene> &genes, bool canonical)
{
    // bed, interval = HALF-OPEN
    //start = 0-based  
    // end   = 1-based  
    // length = end - start

    std::ifstream infile(bedfile);
    if (!infile) {
        std::cerr << "ERROR: could not open file " << bedfile << "\n";
        std::exit(EXIT_FAILURE);
    }
    cout<<"Reading gene annotations from BED file: " << bedfile << "\n";

    std::string line;

    // Temporary storage: gene_name -> vector of entries
    std::map<std::string, std::vector<Gene> > geneMap;

    while (std::getline(infile, line)) {
        if (line.empty()) continue; // skip empty lines

        std::istringstream iss(line);
        std::string chrom, name;
        int start, end;

        if (!(iss >> chrom >> start >> end >> name)) {
            std::cerr << "Warning: malformed/header line skipped -> " << line << "\n";
            continue;
        }

        //check that start is a number, check that start is less than end
        if (start < 0 || end <= start) {
            std::cerr << "Warning: invalid start/end skipped -> " << line << "\n";
            continue;
        }

        // Convert BED 0-based start, 1-based end → internal 1-based inclusive
        int convertedEnd = end;
        int convertedStart = start + 1;


        std::string key = chrom + "|" + name;
        geneMap[key].push_back({chrom, convertedStart, convertedEnd, name});
    }

    infile.close();

    // --- Collapse / deduplicate ---
    genes.clear();
    for (auto &kv : geneMap) {
        auto &entries = kv.second;

        if (canonical) {
            // Pick the longest entry; if multiple longest, pick first
            auto longestIt = std::max_element(entries.begin(), entries.end(),
                                              [](const Gene &a, const Gene &b) {
                                                  return (a.end - a.start) < (b.end - b.start);
                                              });
            genes.push_back(*longestIt);
        } else {
            // Take min start / max end across all entries
            int minStart = entries[0].start;
            int maxEnd   = entries[0].end;
            std::string chrom = entries[0].chrom;

            for (const auto &g : entries) {
                if (g.start < minStart) minStart = g.start;
                if (g.end   > maxEnd)   maxEnd   = g.end;
                // Optional: check all chroms are same
            }
            genes.push_back({chrom, minStart, maxEnd, entries[0].name});
        }
    }
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

double GeneAnalyzer::compute_length_residual(
    const double length,
    const double score,
    double beta0,
    double beta1)
{
    double residual = MISSING_SCORE;

    if (abs(score) != MISSING_SCORE) {
        double pred = beta0 + beta1 * std::log(length);
        residual = score - pred;
    }

    return residual;
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


void GeneAnalyzer::annotateSNPs(std::string geneFile, bool useGTF,
                                  std::vector<std::string> normFiles,
                                  bool XP, string stat,
                                  int minSNPs)
{
            // Force ".norm" prefix behavior (always)
        // {
        //     std::size_t pos = normFile.find(".norm");
        //     if (pos != std::string::npos)
        //     {
        //         normFile = normFile.substr(0, pos + 5); // include ".norm"
        //         std::cout << normFile << "\n";
        //     }
        // }

    cout<<"Reading norm files for annotation...\n";
    
    bool GENERATE_ANNOTATED_NORM = false;

    std::map<std::string, GeneTableEntry> geneTableMap; // aggregate gene scores
    vector<Gene> genes; // all genes from all chromosomes
    if(useGTF)
        readGenesFromGTF(geneFile, genes); // deduplicated genes from gtf, default take gene span
    else
        readGenesFromBed(geneFile, genes); // deduplicated genes from bed, default take gene span

    std::sort(genes.begin(), genes.end(),
              [](const Gene &a, const Gene &b) {
                  if (a.chrom != b.chrom) return a.chrom < b.chrom;
                  if (a.start != b.start) return a.start < b.start;
                  return a.end < b.end;
              });

    // Organize genes by chromosome
    std::unordered_map<std::string, std::vector<GeneNoChr>> genes_by_chr;
    for (const auto &g : genes)
    {
        genes_by_chr[g.chrom].push_back(GeneNoChr{g.start, g.end, g.name});
    }

    genes.clear(); 

    
    for (size_t normFileId = 0; normFileId < normFiles.size(); ++normFileId)
    {
        std::string normFilePath = normFiles[normFileId];
        std::ifstream normFile(normFilePath);
        if (!normFile)
        {
            std::cerr << "ERROR: could not open norm file " << normFilePath << "\n";
            exit(EXIT_FAILURE);
        }

        std::ofstream fout;
        if(GENERATE_ANNOTATED_NORM){
            fout.open(normFilePath + ".ann");

            if (!fout)
            {
                std::cerr << "ERROR: could not open output file " << normFilePath << ".ann\n";
                exit(EXIT_FAILURE);
            }

        }
            

        std::string line;

        std::getline(normFile, line); //header
        if(GENERATE_ANNOTATED_NORM){
            fout << line + "\t" + "overlap_genes" + "\n";
        }
    
        string chr;
        // Parse norm rows:
        // chr id pos freq ihh1 ihh0 ihs norm_ihs crit
        while (std::getline(normFile, line))
        {
            if (line.empty()) continue;
            std::istringstream iss(line);
            std::string chr, id;
            int pos;
            string norm_score_str;
            double norm_score;
            int crit;

            std::vector<std::string> fields;
            std::string token;

            // Tokenize entire line
            while (iss >> token) {
                fields.push_back(token);
            }

            // Basic sanity check
            if (fields.size() < 5) {
                cerr<<"ERROR: invalid line: "<<line<<"\n";
                exit(EXIT_FAILURE);
            }

            // Last two columns (shared semantic meaning)
            try {
                // Mandatory fields (same for all formats)
                chr = fields[0];
                id  = fields[1];
                pos = std::stoi(fields[2]);
                norm_score_str = fields[fields.size() - 2];
                crit       = std::stoi(fields[fields.size() - 1]);
            } catch (const std::exception&) {
                cerr<<"ERROR: invalid line: "<<line<<"\n";
                exit(EXIT_FAILURE);
            }
            
            // if(XP){
            //     //"chr\tcid\tpos\tgpos\tp1\tihh1\tp2\tihh2\txpehh\n"; //+norm_ihs + crit // crit=+1 top crit=-1 bottom
            //     iss >> chr >> junk >> pos >> junk >> junk >> junk >> junk >> junk >> raw_score >> norm_score >> crit;
            // }else if(IHH12){
            //     //"chr\tid\tpos\tp1\tihh12\n";
            //     iss >> chr >> junk >> pos >> junk >> junk >> junk >> crit;
            // }else{
            //     //"chr\tcid\tpos\tfreq\tihh1\tihh0\tihs\n"; //+norm_ihs + crit // crit=+1 top crit=-1 bottom
            //     iss >> chr >> junk >> pos >> junk >> junk >> junk >> raw_score >> norm_score >> crit;
            // }

            // Skip missing norm_ihs
            if (norm_score_str == "NA") continue; // we don't process missing scores
            norm_score = std::stod(norm_score_str);


            //FOR A SINGLE WINDOW (SNP): START [[[
            const std::string &chrom_from_win = chr;

            string overlap_genes = "";
            std::string ov = "";
            auto it = genes_by_chr.find(chrom_from_win);
            if (it == genes_by_chr.end())
            {
                std::cerr << "WARNING: chromosome " << chrom_from_win
                          << " not found in gene BED file.\n";

            }else{
                //chr found in bed file
                const auto &chrGenes = it->second;

                // Find first gene that could overlap position w.start
                size_t gene_idx = 0;
                while (gene_idx < chrGenes.size() && chrGenes[gene_idx].end <= pos)
                    ++gene_idx;

                // For a single-position "window", overlap if gene.start < pos+1 and gene.end > pos
                // i.e., (w.end >= gene.start && w.start < gene.end) with w.start==w.end==pos
                for (size_t j = gene_idx; j < chrGenes.size() && chrGenes[j].start <= pos; ++j)
                {
                    if (pos >= chrGenes[j].start && pos < chrGenes[j].end)
                    {
                        if(GENERATE_ANNOTATED_NORM){
                            if (!ov.empty()) ov += ", ";
                            ov += chrGenes[j].name;
                        }


                        std::string GENE_ID = chrom_from_win + "|" + chrGenes[j].name;
                        auto &gene = geneTableMap[GENE_ID];

                        // first time init
                        if (gene.nWin == 0)
                        {
                            //gene.windowFileId = normFileId;
                            if(XP){
                                //gene.minScore = norm_score;
                                gene.maxScore = norm_score;
                            }else{
                                gene.maxScore = abs(norm_score);
                            }

                            gene.lengthSpan = chrGenes[j].end - chrGenes[j].start + 1; // inclusive span
                            gene.geneEnd = chrGenes[j].end;
                        }
                        else
                        {
                            if(XP){
                                // if (norm_score < gene.minScore)
                                //     gene.minScore = norm_score;
                                
                                if(norm_score > gene.maxScore)
                                    gene.maxScore = norm_score;
                            }else{
                                if (abs(norm_score) > gene.maxScore)
                                    gene.maxScore = abs(norm_score);
                            }
                        }

                        gene.meanScore += norm_score;
                        //gene.scores.push_back(norm_score);
                        gene.nWin += 1;

                        if(crit == 1){
                            gene.n_crit_top += 1.0;
                        }else if(crit == -1){
                            gene.n_crit_bottom += 1.0;
                        }
                    }
                }
            }

            overlap_genes = ov.empty() ? "-" : ov;
            //FOR A SINGLE WINDOW (SNP): END ]]]
            
            // Output annotated line 
            if(GENERATE_ANNOTATED_NORM)
                fout<<line+"\t"+overlap_genes+"\n";
        }

        normFile.close();

        if(GENERATE_ANNOTATED_NORM){
            fout.close();
            std::cout << "Written gene annotated SNPs to " << normFilePath + ".ann\n";
        }
            
    }

    // Gene table printing (kept; assumes fit_length_regression exists)
    {
        std::vector<double> geneLengths;
        std::vector<double> geneScoresMax;
        //std::vector<double> geneScoresMin;


        for (const auto &[gene, row] : geneTableMap)
        {
            if(row.nWin < minSNPs) continue; //ignore genes with too few SNPs

            geneLengths.push_back(row.lengthSpan);
            geneScoresMax.push_back(row.maxScore);
            // if(XP){
            //     geneScoresMin.push_back(row.minScore);
            // }
        }

        std::pair<double,double> beta = fit_length_regression(geneLengths, geneScoresMax);
        
        //std::pair<double,double> beta_min;
        
        // if(XP)
        //     beta_min = fit_length_regression(geneLengths, geneScoresMin);

        std::cout << "Fitted length regression: score = " << beta.first
                  << " + " << beta.second << " * log(length)\n";

        const bool SINGLE_GENE_TABLE = true;
        if (SINGLE_GENE_TABLE)
        {
            std::ofstream genetable;
            std::string genetablefile = geneFile + "." + stat + ".genetable";
            genetable.open(genetablefile.c_str());
            if (genetable.fail())
            {
                std::cerr << "ERROR: " << genetablefile << " " << strerror(errno);
                exit(EXIT_FAILURE);
            }

            genetable << "chr\tstart\tend\tgene\tlen\tnsnps\tmean_score\t"
                         "nsnps_crit\ttop_score\ttop_score_adj\n";

            // if(XP){

            //     genetable << "chr\tstart\tend\tgene\tlen\tnsnps\tmean_score\t"
            //                  "nsnps_crit_top\ttop_score\ttop_score_adj\t"
            //                  "nsnps_crit_bottom\tbottom_score\tbottom_score_adj\n";
            // }else{
                
            // }


            for (const auto &[gene, row] : geneTableMap)
            {
                size_t delim_pos = gene.find("|");
                std::string gene_name = gene.substr(delim_pos + 1);
                std::string chr = gene.substr(0, delim_pos);

                int nwin = row.nWin;          // here: number of SNP records contributing
                double maxScore = row.maxScore;
                double len = row.lengthSpan;
                double mean = row.meanScore / nwin;

                if (nwin < minSNPs) continue;

                double res_score = compute_length_residual(len, maxScore, beta.first, beta.second);
                //double res_score_min = compute_length_residual(len, maxScore, beta_min.first, beta_min.second);
                
                int geneEnd = row.geneEnd;
                int geneStart = geneEnd - static_cast<int>(len) + 1;

                // size_t nScores = row.scores.size();
                // double mean = gsl_stats_mean(row.scores.data(), 1, nScores);
                // double var  = gsl_stats_variance(row.scores.data(), 1, nScores);


                //  if(XP){

                //     genetable << "chr\tstart\tend\tgene\tlen\tnsnps\tmean_score\t"
                //                 "nsnps_crit_top\ttop_score\ttop_score_adj\t"
                //                 "nsnps_crit_bottom\tbottom_score\tbottom_score_adj\n";
                // }else{
                //     genetable << "chr\tstart\tend\tgene\tlen\tnsnps\tmean_score\t"
                //             "nsnps_crit\ttop_score\ttop_score_adj\n";
                // }
                
                genetable << chr << "\t" << geneStart << "\t" << geneEnd << "\t"
                          << gene_name << "\t" << len << "\t" << nwin << "\t"
                          << mean << "\t" << row.n_crit_top << "\t" 
                          << maxScore << "\t" << res_score << "\t" << "\n";

                // if(XP){
                //     genetable << chr << "\t" << geneStart << "\t" << geneEnd << "\t"
                //           << gene_name << "\t" << len << "\t" << nwin << "\t"
                //           << mean << "\t"   
                //           << row.n_crit_top <<"\t" << maxScore << "\t" << res_score << "\t"
                //           << row.n_crit_bottom <<"\t" << row.minScore << "\t" << res_score_min  << "\n";
                // }
                
            }

            genetable.close();
            std::cout << "Written gene table to " << genetablefile << "\n";
        }
    }
}


void GeneAnalyzer::annotateWindows(std::string geneFile, bool useGTF, vector<std::string> windowFiles){
    
    vector<Gene> genes; // all genes from all chromosomes

    if(useGTF)
        readGenesFromGTF(geneFile, genes, true, true); // deduplicated genes from gtf, default take gene span
    else
        readGenesFromBed(geneFile, genes); // deduplicated genes from bed, default take gene span
    //readGenesFromGTF(geneBedFile, genes, true, true); // deduplicated genes from gtf, default take gene span


    std::sort(genes.begin(), genes.end(),
    [](const Gene &a, const Gene &b)
    {
        if (a.chrom != b.chrom)
            return a.chrom < b.chrom; // lexicographic chromosome order
        if (a.start != b.start)
            return a.start < b.start; // sort by start
        return a.end < b.end;         // sort by end
    });

    
    std::unordered_map<std::string, std::vector<GeneNoChr> > genes_by_chr; // Organize genes by chromosome
    for (const auto &g : genes)
    {
        GeneNoChr gnc{g.start, g.end, g.name};
        genes_by_chr[g.chrom].push_back(gnc);
    }

    int numGenesFromAllFiles = genes.size();
    genes.clear();

    std::string windowLine;
    for (size_t windowFileId = 0; windowFileId < windowFiles.size(); ++windowFileId)
    {
        string windowFile = windowFiles[windowFileId];
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

        std::getline(winFile, windowLine); // skip header
        fout << windowLine << "\t" << "overlap_genes"
             << "\n";
        // print header
        // if (XP)
        // {
        //     fout << "chr\tstart\tend\tnSNPs\tfrac_top\tfrac_bottom\tperc_top\tperc_bottom\ttop_score\tbottom_score\toverlap_genes\n";
        // }
        // else
        // {
        //     fout << "chr\tstart\tend\tnSNPs\tfrac_extreme\tperc\tscore\toverlap_genes\n";
        // }


        while (std::getline(winFile, windowLine))
        {
            std::istringstream iss(windowLine);
            //Window w;
            int w_start, w_end;
            std::string w_chrom;

            string overlap_genes="";
            iss >> w_chrom >> w_start >> w_end;
            //insert
            //cout<<"Window: "<< w.chrom << ":" << w.start << "-" << w.end << " Score: " << w.score_max << "\n";
            // check chromosome, only get vector of genes for that chromosome
            
            
            auto it = genes_by_chr.find(w_chrom);
            if (it == genes_by_chr.end())
            {
                std::cerr << "WARNING: chromosome " << w_chrom << " not found in gene BED file.\n";
            }else{
                //chr found in bed file
                const auto &genes = it->second;

                // Two-pointer approach
                size_t gene_idx = 0;
                std::string ov;

                // Advance gene_idx to the first gene that might overlap
                while (gene_idx < genes.size() && genes[gene_idx].end <= w_start)
                {
                    gene_idx++;
                }

                // Check overlapping genes starting from current gene_idx
                size_t j = gene_idx;
                std::string GENE_ID = w_chrom + "|" + genes[j].name;

                while (j < genes.size() && genes[j].start < w_end)
                {
                    if ((w_end > genes[j].start && w_start < genes[j].end))
                    { // overlaps
                        if (!ov.empty())
                            ov += ", ";
                        ov += genes[j].name;
                    }
                    j++;
                }

                overlap_genes = ov == "" ? "-" : ov;
            }
            fout << windowLine << "\t" << overlap_genes << "\n";
        }
        fout.close();
        cout << "Written gene annotated windows to " << windowFile + ".ann\n";
    }
}

// deduplicate genes from GTF/GFF3 by gene_name (or gene_id/Name/ID if not present)
void GeneAnalyzer::readGenesFromGTF(std::string gtffile, std::vector<Gene> &genes,
                                    bool useTranscripts, bool canonical)
{
    igzstream infile;
    std::cerr << "Opening " << gtffile << "...\n";
    infile.open(gtffile.c_str());

    if (!infile) {
        std::cerr << "ERROR: could not open file " << gtffile << "\n";
        std::exit(EXIT_FAILURE);
    }

    auto isGff3Attributes = [](const std::string& attr) -> bool {
        // heuristic: GFF3 looks like key=value;key=value and usually has '='
        // GTF usually has key "value";
        return (attr.find('=') != std::string::npos);
    };

    auto getGtfAttr = [](const std::string& attr, const std::string& key) -> std::string {
        // expects: key "VALUE"
        size_t pos = attr.find(key);
        if (pos == std::string::npos) return "NA";
        size_t firstQuote = attr.find('"', pos);
        if (firstQuote == std::string::npos) return "NA";
        size_t secondQuote = attr.find('"', firstQuote + 1);
        if (secondQuote == std::string::npos) return "NA";
        return attr.substr(firstQuote + 1, secondQuote - firstQuote - 1);
    };

    auto getGff3Attr = [](const std::string& attr, const std::string& key) -> std::string {
        // expects: key=VALUE;key2=VALUE2
        std::string pat = key + "=";
        size_t p = attr.find(pat);
        if (p == std::string::npos) return "NA";
        p += pat.size();
        size_t q = attr.find(';', p);
        std::string val = attr.substr(p, (q == std::string::npos ? attr.size() : q) - p);

        // Some GFF3 values can be URL-encoded; we leave as-is (safe).
        // Also sometimes ID looks like "gene:XYZ"; keep as-is.
        return val.empty() ? "NA" : val;
    };

    auto isProteinCoding = [&](const std::string& attributes) -> bool {
        if (!isGff3Attributes(attributes)) {
            // GTF
            if (attributes.find("gene_biotype \"protein_coding\"") != std::string::npos) return true;
            if (attributes.find("gene_type \"protein_coding\"") != std::string::npos) return true;
            // fallback (looser)
            if (attributes.find("\"protein_coding\"") != std::string::npos) return true;
            return false;
        } else {
            // GFF3
            if (attributes.find("biotype=protein_coding") != std::string::npos) return true;
            if (attributes.find("gene_biotype=protein_coding") != std::string::npos) return true;
            if (attributes.find("gene_type=protein_coding") != std::string::npos) return true;
            return false;
        }
    };

    std::string line;
    // key -> vector of spans {chrom, start, end, geneName}
    std::map<std::string, std::vector<Gene>> geneMap;

    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string chrom, source, feature, score, strand, frame;
        int start, end;

        if (!(iss >> chrom >> source >> feature >> start >> end >> score >> strand >> frame))
            continue;

        std::string attributes;
        std::getline(iss, attributes); // remainder (starts with whitespace)

        // --- Decide format based on attributes ---
        const bool isGFF3 = isGff3Attributes(attributes);

        // --- Filter feature type ---
        if (!useTranscripts) {
            if (feature != "gene") continue;
            // Only keep protein-coding genes
            if (!isProteinCoding(attributes)) continue;
        } else {
            // GTF: transcript, GFF3: mRNA is common (also "transcript" sometimes)
            if (!(feature == "transcript" || feature == "mRNA")) continue;
            // Optional: you *may* still want protein-coding only here too:
            // if (!isProteinCoding(attributes)) continue;
        }

        // --- Extract gene identifier/name ---
        std::string geneName = "NA";
        if (!isGFF3) {
            // GTF
            geneName = getGtfAttr(attributes, "gene_name");
            if (geneName == "NA") geneName = getGtfAttr(attributes, "gene_id");
        } else {
            // GFF3
            geneName = getGff3Attr(attributes, "Name");
            if (geneName == "NA") geneName = getGff3Attr(attributes, "gene_name");
            if (geneName == "NA") {
                // fall back to ID
                geneName = getGff3Attr(attributes, "ID");
            }
        }

        if (geneName == "NA") continue;

        // Convert coordinates:
        // GTF and GFF3 are 1-based inclusive.
        // Your internal representation appears to be 1-based inclusive too, so keep as-is.
        int convertedStart = start;
        int convertedEnd   = end;

        std::string key = chrom + "|" + geneName;
        geneMap[key].push_back({chrom, convertedStart, convertedEnd, geneName});
    }

    infile.close();

    // --- Collapse / deduplicate ---
    genes.clear();
    genes.reserve(geneMap.size());

    for (auto &kv : geneMap) {
        auto &entries = kv.second;
        if (entries.empty()) continue;

        if (canonical) {
            // Pick the longest entry; tie-breaker = first
            auto longestIt = std::max_element(entries.begin(), entries.end(),
                                              [](const Gene &a, const Gene &b) {
                                                  return (a.end - a.start) < (b.end - b.start);
                                              });
            genes.push_back(*longestIt);
        } else {
            // Take min start / max end across all entries
            int minStart = entries[0].start;
            int maxEnd   = entries[0].end;
            std::string chrom0 = entries[0].chrom;
            std::string geneName0 = entries[0].name;

            for (const auto &g : entries) {
                if (g.start < minStart) minStart = g.start;
                if (g.end   > maxEnd)   maxEnd   = g.end;
            }
            genes.push_back({chrom0, minStart, maxEnd, geneName0});
        }
    }
}



// // deduplicate genes from GTF by gene_name (or gene_id if gene_name not present)
// void GeneAnalyzer::readGenesFromGTF(std::string gtffile, std::vector<Gene> &genes,
//                                     bool useTranscripts, bool canonical)
// {
//     igzstream infile;
//     std::cerr << "Opening " << gtffile << "...\n";
//     infile.open(gtffile.c_str());

//     if (!infile) {
//         std::cerr << "ERROR: could not open file " << gtffile << "\n";
//         std::exit(EXIT_FAILURE);
//     }

//     std::string line;

//     // Temporary storage: gene_name -> vector of spans {chrom, start, end}
//     std::map<std::string, std::vector<Gene>> geneMap;

//     while (std::getline(infile, line)) {
//         if (line.empty() || line[0] == '#') continue;

//         std::istringstream iss(line);
//         std::string chrom, source, feature, score, strand, frame, attributes;
//         int start, end;

//         if (!(iss >> chrom >> source >> feature >> start >> end >> score >> strand >> frame))
//             continue;

//         std::getline(iss, attributes);

//         // Extract gene_name
//         std::string geneName = "NA";
//         size_t pos = attributes.find("gene_name");
//         if (pos != std::string::npos) {
//             size_t firstQuote = attributes.find('"', pos);
//             size_t secondQuote = attributes.find('"', firstQuote + 1);
//             geneName = attributes.substr(firstQuote + 1, secondQuote - firstQuote - 1);
//         } else {
//             pos = attributes.find("gene_id");
//             if (pos != std::string::npos) {
//                 size_t firstQuote = attributes.find('"', pos);
//                 size_t secondQuote = attributes.find('"', firstQuote + 1);
//                 geneName = attributes.substr(firstQuote + 1, secondQuote - firstQuote - 1);
//             }
//         }

//         // Convert GTF (1-based inclusive) //internal 1 based inclusive same // internal (0-based inclusive)
//         int convertedStart = start; // - 1;
//         int convertedEnd   = end; // - 1;

//         if (!useTranscripts) {
//             if (!(attributes.find("\"protein_coding\"") != std::string::npos)) {
//                 continue;
//             }
//             if (feature != "gene") continue;
//         } else {
//             if (feature != "transcript") continue;
//         }

//         string key = chrom + "|" + geneName;
//         geneMap[key].push_back({chrom, convertedStart, convertedEnd, geneName});
//     }

//     infile.close();

//     // --- Collapse / deduplicate ---
//     genes.clear();
//     for (auto &kv : geneMap) {
//         auto &entries = kv.second;

//         if (canonical) {
//             // Pick the **longest entry**; if multiple longest, pick the first
//             auto longestIt = std::max_element(entries.begin(), entries.end(),
//                                               [](const Gene &a, const Gene &b) {
//                                                   return (a.end - a.start) < (b.end - b.start);
//                                               });
//             genes.push_back(*longestIt);
//         } else {
//             // Take **min start / max end** across all entries
//             int minStart = entries[0].start;
//             int maxEnd   = entries[0].end;
//             std::string chrom = entries[0].chrom;

//             for (const auto &g : entries) {
//                 if (g.start < minStart) minStart = g.start;
//                 if (g.end   > maxEnd)   maxEnd   = g.end;
//                 // Optional: check all chroms are same
//             }
//             genes.push_back({chrom, minStart, maxEnd, kv.first});
//         }
//     }
// }
