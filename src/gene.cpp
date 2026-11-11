
#include "gene.h"
#include <unordered_map>
#include "gzstream.h"

// deduplicated genes from bed, default take gene span
void GeneAnalyzer::readGenesFromBed(std::string bedfile, std::vector<Gene> &genes, bool canonical)
{
    // bed, interval = HALF-OPEN
    // start = 0-based  
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

void GeneAnalyzer::perm_test(string fileA, string fileB, long seed) {
    // norm's gene tables carry the length-adjusted maximum as top_score_adj
    // (3.0.1 looked for a column "max_lenreg", which no table has; B38)
    std::vector<double> scoresA = read_scores(fileA, "gene", "top_score_adj");
    std::vector<double> scoresB = read_scores(fileB, "gene", "top_score_adj");
    if (scoresA.empty() || scoresB.empty()) {
        std::cerr << "ERROR: the permutation test needs at least one gene with a numeric top_score_adj in each set ("
                  << scoresA.size() << " in " << fileA << ", " << scoresB.size() << " in " << fileB << ").\n";
        exit(EXIT_FAILURE);
    }

    // --- Permutation test (one-sided: is the target set's mean higher?) ---
    int n_perms = 10000;
    int count = 0;
    std::vector<double> combined = scoresA;
    combined.insert(combined.end(), scoresB.begin(), scoresB.end());

    int nA = scoresA.size();
    int n_total = combined.size();
    std::vector<int> indices(n_total);
    std::iota(indices.begin(), indices.end(), 0);

    if (seed < 0) seed = static_cast<long>(std::time(nullptr));
    snum::Rng rng(static_cast<uint64_t>(seed));   // was gsl_rng_default, so a seed gives different permutations than GSL builds

    double T_obs = snum::mean(scoresA.data(), scoresA.size()) -
                snum::mean(scoresB.data(), scoresB.size());

    for (int i = 0; i < n_perms; ++i) {
        rng.shuffle(indices.data(), n_total);

        std::vector<double> permA(nA), permB(n_total - nA);
        for (int j = 0; j < nA; ++j) permA[j] = combined[indices[j]];
        for (int j = nA; j < n_total; ++j) permB[j - nA] = combined[indices[j]];

        double T_perm = snum::mean(permA.data(), permA.size()) -
                        snum::mean(permB.data(), permB.size());
        if (T_perm >= T_obs) count++;
    }

    // add-one estimate: the observed labelling is one of the permutations,
    // so the p-value is never 0
    double p_value = static_cast<double>(count + 1) / (n_perms + 1);
    std::cout << "Target genes: " << nA << ", background genes: " << (n_total - nA) << "\n";
    std::cout << "Permutations: " << n_perms << ", seed: " << seed << " (--perm-seed)\n";
    std::cout << "Observed mean difference (A - B): " << T_obs << "\n";
    std::cout << "Empirical one-sided p-value (A > B): " << p_value << "\n";
    std::cout << ((p_value < 0.05) ? "Set A is significantly higher than Set B.\n"
                                : "No significant difference between Set A and Set B.\n");

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

std::pair<double,double> GeneAnalyzer::fit_length_regression(
    const std::vector<double> &lengths,
    const std::vector<double> &scores)
{
    // score ~ intercept + log(length), over genes with a score
    std::vector<double> x, y;
    x.reserve(lengths.size());
    y.reserve(lengths.size());
    for (size_t i = 0; i < lengths.size(); i++) {
        if (abs(scores[i]) != MISSING_SCORE) {
            x.push_back(std::log(lengths[i]));
            y.push_back(scores[i]);
        }
    }
    // With 2 genes every residual is 0. Callers require at least 3 (B25).
    if (x.size() <= 2) return {0.0, 0.0};
    return snum::fit_line(x, y);
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
        struct SnpScan { size_t idx = 0; long lastPos = -1; };
        std::unordered_map<std::string, SnpScan> snpScan;   // P3: per-chromosome gene scan position
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

                // Find first gene that could overlap position w.start.
                // P3: genes are sorted by start and SNPs by position, and a
                // gene ending before this SNP cannot overlap a later one, so
                // the scan resumes where the previous SNP's scan stopped
                // (reset if positions go backwards). 3.0.1 restarted at gene
                // 0 for every SNP.
                SnpScan &sc = snpScan[chrom_from_win];
                if (pos < sc.lastPos) sc.idx = 0;
                sc.lastPos = pos;
                size_t gene_idx = sc.idx;
                // genes are 1-based inclusive [start, end] (B37)
                while (gene_idx < chrGenes.size() && chrGenes[gene_idx].end < pos)
                    ++gene_idx;
                sc.idx = gene_idx;

                // For a single-position "window", overlap if gene.start < pos+1 and gene.end > pos
                // i.e., (w.end >= gene.start && w.start < gene.end) with w.start==w.end==pos
                for (size_t j = gene_idx; j < chrGenes.size() && chrGenes[j].start <= pos; ++j)
                {
                    if (pos >= chrGenes[j].start && pos <= chrGenes[j].end)
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
                                gene.minScore = norm_score;
                                gene.maxScore = norm_score;
                            }else{
                                double score_to_compare = (stat=="ihh12") ? norm_score : abs(norm_score);
                                gene.maxScore = score_to_compare;
                            }

                            gene.lengthSpan = chrGenes[j].end - chrGenes[j].start + 1; // inclusive span
                            gene.geneEnd = chrGenes[j].end;
                        }
                        else
                        {
                            if(XP){
                                if (norm_score < gene.minScore)
                                    gene.minScore = norm_score;
                                
                                if(norm_score > gene.maxScore)
                                    gene.maxScore = norm_score;
                            }else{
                                double score_to_compare = (stat=="ihh12") ? norm_score : abs(norm_score);
                                if(score_to_compare > gene.maxScore)
                                    gene.maxScore = score_to_compare;
                            }
                        }

                        //gene.meanScore += norm_score;
                        //gene.scores.push_back(norm_score);
                        gene.nWin += 1;

                        if(crit == 1){
                            gene.nsnps_crit_top += 1.0;
                        }else if(crit == -1){
                            gene.nsnps_crit_bot += 1.0;
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
        std::vector<double> geneScoresMin;


        for (const auto &[gene, row] : geneTableMap)
        {
            if(row.nWin < minSNPs) continue; //ignore genes with too few SNPs

            geneLengths.push_back(row.lengthSpan);
            geneScoresMax.push_back(row.maxScore);
            if(XP){
                geneScoresMin.push_back(row.minScore);
            }
        }

        // The length regression has two coefficients; with fewer than three
        // genes it is undetermined or fits exactly, and GSL aborts below two
        // (B25). Skip it and write NA in the *_adj columns.
        const bool FIT_LENGTH = geneLengths.size() >= 3;
        std::pair<double,double> beta = {0.0, 0.0};
        std::pair<double,double> beta_min = {0.0, 0.0};

        if (FIT_LENGTH) {
            beta = fit_length_regression(geneLengths, geneScoresMax);
            if(XP)
                beta_min = fit_length_regression(geneLengths, geneScoresMin);
            std::cout << "Fitted length regression: score = " << beta.first
                      << " + " << beta.second << " * log(length)\n";
        } else {
            std::cerr << "WARNING: " << geneLengths.size() << " gene(s) with at least " << minSNPs
                      << " SNPs; at least 3 are needed for the length regression. "
                      << "It is skipped and the *_adj columns are NA.\n";
        }

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



            if(XP){

                genetable << "chr\tstart\tend\tgene\tlen\tnsnps\t"
                             "frac_crit_top\ttop_score\ttop_score_adj\t"
                             "frac_crit_bottom\tbottom_score\tbottom_score_adj\n";
            }else{
                genetable << "chr\tstart\tend\tgene\tlen\tnsnps\t"
                         "frac_crit\ttop_score\ttop_score_adj\n"; 
            }  


            for (const auto &[gene, row] : geneTableMap)
            {
                size_t delim_pos = gene.find("|");
                std::string gene_name = gene.substr(delim_pos + 1);
                std::string chr = gene.substr(0, delim_pos);

                int nwin = row.nWin;          // here: number of SNP records contributing
                double maxScore = row.maxScore;
                double len = row.lengthSpan;
                
                //double mean = row.meanScore / nwin;

                if (nwin < minSNPs) continue;

                std::string res_score = "NA", res_score_min = "NA";
                if (FIT_LENGTH) {
                    std::ostringstream a, b;
                    a << compute_length_residual(len, maxScore, beta.first, beta.second);
                    // bottom score against the bottom-score fit (was maxScore, B26)
                    b << compute_length_residual(len, row.minScore, beta_min.first, beta_min.second);
                    res_score = a.str(); res_score_min = b.str();
                }
                
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
                
                

                if(XP){
                    genetable << chr << "\t" << geneStart << "\t" << geneEnd << "\t"
                          << gene_name << "\t" << len << "\t" << nwin << "\t"
                          << row.nsnps_crit_top*1.0/nwin <<"\t" << maxScore << "\t" << res_score << "\t"
                          << row.nsnps_crit_bot*1.0/nwin <<"\t" << row.minScore << "\t" << res_score_min  << "\n";
                }else{
                    genetable << chr << "\t" << geneStart << "\t" << geneEnd << "\t"
                          << gene_name << "\t" << len << "\t" << nwin << "\t"
                          << row.nsnps_crit_top*1.0/nwin << "\t" 
                          << maxScore << "\t" << res_score << "\t" << "\n";

                }
                
            }

            genetable.close();
            std::cout << "Written gene table to " << genetablefile << "\n";
        }
    }
}


void GeneAnalyzer::annotateWindows(std::string geneFile, bool useGTF, vector<std::string> windowFiles){
    
    vector<Gene> genes; // all genes from all chromosomes

    if(useGTF)
        readGenesFromGTF(geneFile, genes); // deduplicated genes from gtf, default take gene span
    else
        readGenesFromBed(geneFile, genes); // deduplicated genes from bed, default take gene span
    //readGenesFromGTF(geneBedFile, genes, true, true); // deduplicated genes from gtf, default take gene span


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
    // int numGenesFromAllFiles = genes.size();
    // genes.clear();

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


        struct WinScan { size_t idx = 0; long lastStart = -1; };
        std::unordered_map<std::string, WinScan> winScan;   // P3: per-chromosome gene scan position
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
                const auto &chrGenes = it->second;

                // Two-pointer approach (P3: resumes where the previous
                // window's scan stopped; reset if windows go backwards)
                WinScan &ws = winScan[w_chrom];
                if (w_start < ws.lastStart) ws.idx = 0;
                ws.lastStart = w_start;
                size_t gene_idx = ws.idx;
                std::string ov;

                // Advance gene_idx to the first gene that might overlap
                // genes and windows are both 1-based inclusive (B37)
                while (gene_idx < chrGenes.size() && chrGenes[gene_idx].end < w_start)
                {
                    gene_idx++;
                }
                ws.idx = gene_idx;

                // Check overlapping genes starting from current gene_idx
                size_t j = gene_idx;
                //std::string GENE_ID = w_chrom + "|" + chrGenes[j].name;

                while (j < chrGenes.size() && chrGenes[j].start <= w_end)
                {
                    if (w_end >= chrGenes[j].start && w_start <= chrGenes[j].end)
                    { // overlaps
                        if (!ov.empty())
                            ov += ", ";
                        ov += chrGenes[j].name;
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

// deduplicate genes from GTF/GFF3 by gene_name (or gene_id/Name/ID if not present) - only protein-coding genes
void GeneAnalyzer::readGenesFromGTF(std::string gtffile, std::vector<Gene> &genes,
                                    bool useTranscripts, bool canonical)
{
    igzstream infile;
    infile.open(gtffile.c_str());

    if (!infile) {
        std::cerr << "ERROR: could not open file " << gtffile << "\n";
        std::exit(EXIT_FAILURE);
    }
    std::cout<<"Reading gene annotations from GTF file: " << gtffile << "\n";

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
