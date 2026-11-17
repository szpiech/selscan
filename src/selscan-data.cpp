
#include "selscan-data.h"
#include "gzlinereader.h"
#include <set>

void HapMap::initParamsInHap(HapData &hapData){
    hapData.MISSING_ALLOWED = p.MISSING_ALLOWED;
    hapData.unphased = p.UNPHASED;
    hapData.MULTI_MAF = p.MULTI_MAF;

    //hapData.MULTI_CHR = p.MULTI_CHR;
    // if(p.MULTI_CHR){
    //     cerr<<"Multi-chromosome support enabled."<<endl;
    //     // MULTI-chromosome support added
        
    //     if(p.CHR_LIST!=""){
    //         std::stringstream ss(p.CHR_LIST);
    //         std::string item;
    //         while (std::getline(ss, item, ',')) { // Split the input string by commas
    //             if(item!="")
    //                 chr_set.insert(item);
    //         }
    //     }

    //     if(chr_set.empty()){
    //         cerr<<"ERROR: No valid chromosome list provided.\n";
    //     }
    // }

    // if we are in multi param mode, or using xp, we want to set the MAF cutoff to the minimum
    this->MIN_MAF_CUTOFF = p.MAF;

    if(p.CALC_IHS || p.CALC_NSL || p.CALC_SOFT){
        if(p.MULTI_MAF){
            for (const auto& param : ps) {
                if(param.MAF < MIN_MAF_CUTOFF){
                    this->MIN_MAF_CUTOFF = param.MAF;
                }
            }
        }
        // A co-requested statistic that uses every site (pi, EHH, EHH12, XP)
        // must see the unfiltered data, so the loader cannot prefilter. The
        // MAF filter for iHS/nSL/iHH12 is then applied at run time by their
        // kernels (MULTI_MAF is set for every combined run), which reproduces
        // the single-run result. Without this, pi and EHH in a combined run
        // were computed on MAF-filtered sites (B17).
        if(p.CALC_PI || p.SINGLE_EHH || p.SINGLE_EHH12 || p.CALC_XP || p.CALC_XPNSL){
            this->MIN_MAF_CUTOFF = 0;
        }
    }else{
        // if we are not calculating IHS or NSL, we can set the MAF cutoff to 0
        this->MIN_MAF_CUTOFF = 0;
    }
    // if IHS OR NSL othwerwise MINMAF==0
}




// ms hap format : number of rows = number of haplotypes, number of columns = number of loci
void HapMap::readHapDataMSHap(string filename, HapData &hapData)
{
    //HANDLE_ERROR("Bug found in ms-style hap format. This function is undergoing development.");
    initParamsInHap(hapData);

    // Phase === PHASE 1: Determine nloci, nhaps, and MAF-based skiplist ===
    //PHASE 1: Read input file to get "nloci", "nhaps" and "skiplist"
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        exit(EXIT_FAILURE);
    }

    // Each row is one haplotype. With --unphased, rows 2i and 2i+1 are the
    // two haplotypes of individual i. 3.0.1 handled rows in pairs in both
    // modes, so with an odd row count the last haplotype was never read or
    // checked and stayed all-ancestral (B14). Blank lines are skipped (U9).
    auto stripSpace = [](string &s) {
        s.erase(std::remove_if(s.begin(), s.end(), [](char c) {
            return std::isspace(static_cast<unsigned char>(c));
        }), s.end());
    };

    string line;
    string prev_line;
    int current_nhaps = 0;

    vector<int> skiplist;

    int lineNum = 0;
    while (getline(fin, prev_line)) { lineNum++; stripSpace(prev_line); if (!prev_line.empty()) break; }
    if (prev_line.empty()) {
        HANDLE_ERROR(filename + " has no haplotype rows.");
    }

    int nloci_before_filter = prev_line.length(); // number of loci before MAF filtering

    vector<int> number_of_1s_per_locus(nloci_before_filter, 0);
    vector<int> number_of_2s_per_locus(nloci_before_filter, 0);

    // PHASE 1: count rows and alleles, and check the format, before anything
    // is allocated.
    int num_rows = 0;
    line = prev_line;
    do
    {
        stripSpace(line);
        if (line.empty()) continue;

        if (static_cast<int>(line.length()) != nloci_before_filter) {
            HANDLE_ERROR("line " + std::to_string(lineNum) + " of " + filename + " has " + std::to_string(line.length()) +
                          " loci, but the first line has " + std::to_string(nloci_before_filter) + ".");
        }
        for (int loc = 0; loc < nloci_before_filter; loc++) {
            if (line[loc] != '0' && line[loc] != '1') {
                HANDLE_ERROR("line " + std::to_string(lineNum) + " of " + filename + " has '" + string(1, line[loc]) +
                              "' at locus " + std::to_string(loc + 1) + ". Alleles must be coded 0/1 only.");
            }
        }

        if (p.UNPHASED) {
            if (num_rows % 2 == 1) {
                for (int loc = 0; loc < nloci_before_filter; loc++) {
                    int g = (line[loc] == '1') + (prev_line[loc] == '1');
                    if (g == 2) number_of_2s_per_locus[loc]++;
                    else if (g == 1) number_of_1s_per_locus[loc]++;
                }
            }
        } else {
            for (int loc = 0; loc < nloci_before_filter; loc++) {
                if (line[loc] == '1') number_of_1s_per_locus[loc]++;
            }
        }

        prev_line = line;
        num_rows++;
    } while ((lineNum++, getline(fin, line)));

    if (p.UNPHASED && num_rows % 2 != 0) {
        HANDLE_ERROR(filename + " has " + std::to_string(num_rows) +
                      " haplotype rows. With " + ARG_UNPHASED + " the rows are read as pairs (two per individual), so the number must be even.");
    }
    current_nhaps = (p.UNPHASED) ? num_rows/2 : num_rows;

    int n_maf_skipped = 0;
    for(int loc=0; loc<nloci_before_filter; loc++){
        if(preSkipAt(loc)){
            skiplist.push_back(loc);
        }else if(shouldSkipLocus(number_of_1s_per_locus[loc], number_of_2s_per_locus[loc], num_rows)){
            skiplist.push_back(loc);
            n_maf_skipped++;
        }
    }

    fin.clear(); // clear error flags
    fin.close();

    number_of_1s_per_locus.clear();
    number_of_2s_per_locus.clear();

    // === Logging summary ===

    if (p.SKIP) {
        LOG(ARG_SKIP << " set. Removing all variants < " << p.MAF);
        LOG("Removed " << n_maf_skipped << " low frequency variants from haplotype data.");
    } else {
        LOG(ARG_KEEP << " set. NOT removing variants < " << p.MAF);
    }

    hapData.initHapData(current_nhaps, nloci_before_filter - skiplist.size());

    //PHASE 2: Open input File 2nd time To Load into Data Structure
    LOG("Re-loading input to read " << current_nhaps << " haplotypes and "
        << (nloci_before_filter - skiplist.size()) << " loci into data structures.");

    fin.open(filename.c_str());
    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        exit(EXIT_FAILURE);
    }

    prev_line = "";
    int hap = 0;
    while (hap < num_rows && getline(fin, line))
    {
        stripSpace(line);
        if (line.empty()) continue;

        int skip_idx = 0;
        int loc_after_filter = 0;
        for(int loc=0; loc<nloci_before_filter; loc++){

            if (skip_idx < static_cast<int>(skiplist.size()) && loc == skiplist[skip_idx]) {
                if(hap==num_rows-1){
                    hapData.skipQueue.push(loc);
                }
                ++skip_idx;  // move to the next skip index
                continue;
            }

            if(p.UNPHASED){
                if (hap % 2 == 1){
                    int g = (line[loc] == '1') + (prev_line[loc] == '1');
                    if (g == 2) hapData.addAllele2(loc_after_filter, (hap-1)/2);
                    else if (g == 1) hapData.addAllele1(loc_after_filter, (hap-1)/2);
                }
            }else{
                if (line[loc]=='1') hapData.addAllele1(loc_after_filter, hap);
            }
            loc_after_filter++;
        }

        prev_line = line; // store the previous line
        hap++;
    }
    fin.close();

}




// B33: shape of a GT field. Returns 0 for a diploid call ('a|b' or 'a/b'
// followed by ':', a tab or the end of the line), 1 for a haploid call ('a'
// alone) and 2 for anything else. 3.0.1 read a fixed 3 characters, so a
// haploid call took its second allele from the next sample's field.
static inline int gtShape(const char* g, const char* end)
{
    const ptrdiff_t n = end - g;
    // common case first: 'a|b' or 'a/b' then a tab, ':' or the end of the line
    if (n >= 3 && (g[1] == '|' || g[1] == '/') && g[2] != '\t' && g[2] != ':' &&
        (n == 3 || g[3] == '\t' || g[3] == ':' || g[3] == '\r' || g[3] == '\n'))
        return (g[0] == '\t' || g[0] == ':') ? 2 : 0;
    auto stop = [end](const char* q) { return q >= end || *q == '\t' || *q == ':' || *q == '\r' || *q == '\n'; };
    if (stop(g)) return 2;
    if (stop(g + 1)) return 1;
    return 2;
}

// text of the genotype field at g, for error messages
static inline string gtText(const char* g, const char* end)
{
    const char* q = g;
    while (q < end && *q != '\t' && *q != ':' && *q != '\r' && *q != '\n') ++q;
    return string(g, q);
}

// true if the line has a non-whitespace character (U9: blank lines are skipped)
static inline bool hasData(const string &l)
{
    for (char c : l) if (!std::isspace(static_cast<unsigned char>(c))) return true;
    return false;
}

// 0/1 coding of one row of a --thap or --tped file, checked in pass 1 for
// every row and in both modes (3.0.1 never checked rows with --unphased).
// With --unphased the haplotype count is checked on the first row, before
// the pairs are read.
void HapMap::checkHapRow(const string &haps, const string &filename, int lineNum, bool firstRow)
{
    for (size_t k = 0; k < haps.size(); k++) {
        if (haps[k] != '0' && haps[k] != '1') {
            HANDLE_ERROR("line " + std::to_string(lineNum) + " of " + filename + " has '" + string(1, haps[k]) +
                          "' for haplotype " + std::to_string(k + 1) + ". Alleles must be coded 0/1 only.");
        }
    }
    if (firstRow && p.UNPHASED && haps.size() % 2 != 0) {
        HANDLE_ERROR(filename + " has " + std::to_string(haps.size()) + " haplotypes. With " + ARG_UNPHASED +
                      " they are read as pairs (two per individual), so the number must be even.");
    }
}

//reads in haplotype data and also does basic checks on integrity of format
//returns a populated HaplotypeData structure if successful
//impute hap IMPUTE HAP is transposed format where row represents loci,  column replesent individual
//so wc -l of impute hap is same as wc -l of map.
void HapMap::readHapDataTHAP(string filename, HapData& hapData)
{
    initParamsInHap(hapData);

    //PHASE 1: Read input file to get "nloci", "nhaps" and "skiplist"
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        HANDLE_ERROR("Failed to open " + filename + " for reading.");
    }

    string line;
    
    int previous_nhaps = -1;
    int current_nhaps = 0;

    queue<int> skiplist;
    int nloci_before_filter = 0;
    int n_maf_skipped = 0;

    //Counts number of haps (rows) and number of loci (cols)
    //if any lines differ, send an error message and throw an exception
    int lineNum = 0;
    while (getline(fin, line))
    {
        lineNum++;
        if (!hasData(line)) continue;   // U9
        nloci_before_filter++;
        line.erase(std::remove_if(line.begin(), line.end(), [](char c) {
            return std::isspace(static_cast<unsigned char>(c));
        }), line.end());
        current_nhaps = line.length(); // number of haps 
        checkHapRow(line, filename, lineNum, previous_nhaps < 0);

        int number_of_2s = 0; // for unphased
        int number_of_1s = 0;
        for (int hap = 0; hap < current_nhaps; hap++){ 
            if(p.UNPHASED){
                if (hap % 2 == 0){
                    if(line[hap]=='1'  && line[hap+1]=='1'){
                        number_of_2s++;
                    }else if ( (line[hap]=='1' && line[hap+1]=='0') || (line[hap]=='0' && line[hap+1]=='1') ){ // ==1
                        number_of_1s++;
                    }
                }
            }else{
                if(line[hap]=='1'){
                    number_of_1s++;
                }
            }
        }

        if(preSkipAt(nloci_before_filter-1)){
            skiplist.push(nloci_before_filter-1);
        }else if(shouldSkipLocus(number_of_1s, number_of_2s, current_nhaps)){
            skiplist.push(nloci_before_filter-1);
            n_maf_skipped++;
        }
        
    
        if (previous_nhaps < 0)
        {
            previous_nhaps = current_nhaps;
            continue;
        }
        else if (previous_nhaps != current_nhaps) {
            HANDLE_ERROR(("line " + std::to_string(nloci_before_filter) +
                          " of " + filename + " has " + std::to_string(current_nhaps) +
                          " haps, but the previous line has " + std::to_string(previous_nhaps) + ".").c_str());
        }
        
        previous_nhaps = current_nhaps;
    }
    fin.clear(); // clear error flags
    fin.close();

    if(p.SKIP) { //prefilter all sites < MAF
        LOG(ARG_SKIP << " set. Removing all variants < " << p.MAF << ".");
        LOG("Removed "<<n_maf_skipped <<" low frequency variants.");

    }else{
        LOG(ARG_KEEP << " set. NOT removing variants < " << p.MAF << ".");
    }
    
    if(p.UNPHASED){
        if (current_nhaps % 2 != 0)
        {
            HANDLE_ERROR("Number of haplotypes must be even for unphased.");
        }
    }
    current_nhaps = (p.UNPHASED) ? current_nhaps/2 : current_nhaps;
    
    LOG("Loading " << current_nhaps << " haplotypes and " << nloci_before_filter - skiplist.size() << " loci...");
    hapData.initHapData(current_nhaps, nloci_before_filter-skiplist.size());
    hapData.skipQueue = queue<size_t>();

    //PHASE 2: Open input File 2nd time To Load into Data Structure
    fin.open(filename.c_str());
    if (fin.fail())
    {
        HANDLE_ERROR("Failed to open " + filename + " for reading.");
    }
    int locus_after_filter = 0;
    for (int locus = 0; locus < nloci_before_filter; locus++)
    {
        while (getline(fin, line) && !hasData(line)) {}   // U9
        if(!skiplist.empty()){
            if(skiplist.front() == locus){
                skiplist.pop();
                hapData.skipQueue.push(locus);
                continue;
            }
        }

        line.erase(std::remove_if(line.begin(), line.end(), [](char c) {
                return std::isspace(static_cast<unsigned char>(c));
            }), line.end());

        int cols = (p.UNPHASED) ? current_nhaps*2 : current_nhaps;

        if (p.UNPHASED){
            
            
            for (int hap = 0; hap < cols; hap++){ 
                if (hap % 2 == 0){
                    if(line[hap]=='1'  && line[hap+1]=='1'){
                        hapData.addAllele2(locus_after_filter, hap/2);
                    }else if ( (line[hap]=='1' && line[hap+1]=='0') || (line[hap]=='0' && line[hap+1]=='1') ){ // ==1
                        hapData.addAllele1(locus_after_filter, hap/2);
                    }
                }
            }
        }else{ // PHASED (0/1 coding checked in pass 1)
            for (int hap = 0; hap < current_nhaps; hap++)
            {
                if(line[hap]=='1'){
                    hapData.addAllele1(locus_after_filter, hap);
                }
            }
        }

        locus_after_filter++;
    }
    fin.close();
}



// Column 4 (physical position) of every row of a TPED or map file.
static vector<long> readPositionColumn(const string& filename)
{
    GzLineReader fin;   // P1: buffered reader
    fin.open(filename.c_str());
    if (fin.fail()) {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 1;
    }
    vector<long> pos;
    string line, a, b, c, f4;
    int lineNum = 0;
    while (getline(fin, line)) {
        lineNum++;
        if (!hasData(line)) continue;   // U9
        istringstream ss(line);
        int d = 0;
        if (!(ss >> a >> b >> c >> f4)) {
            cerr << "ERROR: line " << lineNum << " of " << filename << " has no physical position in column 4.\n";
            throw 1;
        }
        if (!parsePosition(f4, d)) {   // B44
            cerr << "ERROR: position '" << f4 << "' on line " << lineNum << " of " << filename << POS_LIMIT_MSG << "\n";
            throw 1;
        }
        pos.push_back(d);
    }
    fin.close();
    return pos;
}

// Derived-allele count per locus, and alleles per locus, for TPED ('t'),
// thap ('h': one line per locus) or ms-style hap ('m': one line per
// haplotype) input. Every allele field or character is one allele; '1' is
// derived.
static void scanDerivedCounts(const string& filename, char fmt, vector<long>& derived, long& nalleles)
{
    GzLineReader fin;   // P1: buffered reader
    fin.open(filename.c_str());
    if (fin.fail()) {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 1;
    }
    derived.clear(); nalleles = 0;
    string line, tok;
    while (getline(fin, line)) {
        if (!hasData(line)) continue;   // U9
        if (fmt == 't') {
            istringstream ss(line);
            long n = 0, d = 0; int col = 0;
            while (ss >> tok) { if (++col <= 4) continue; n++; d += (tok == "1"); }
            derived.push_back(d); nalleles = n;
        } else {
            line.erase(std::remove_if(line.begin(), line.end(), [](char c) { return std::isspace(static_cast<unsigned char>(c)); }), line.end());
            if (fmt == 'h') {
                derived.push_back(std::count(line.begin(), line.end(), '1')); nalleles = line.size();
            } else {
                if (derived.empty()) derived.assign(line.size(), 0);
                for (size_t i = 0; i < line.size() && i < derived.size(); i++) derived[i] += (line[i] == '1');
                nalleles++;
            }
        }
    }
    fin.close();
}

// Drop every record at a position that occurs more than once (PI decision):
// a duplicated position cannot be placed on the map, and keeping one of the
// records would be arbitrary. This is what the VCF readers already did; the
// TPED and hap/thap readers kept duplicates.
void HapMap::markDuplicatePositions(const vector<long>& positions, const string& source)
{
    std::map<long, int> count;
    for (long x : positions) count[x]++;
    preSkip.assign(positions.size(), 0);
    size_t nrec = 0, npos = 0;
    for (const auto& kv : count) if (kv.second > 1) npos++;
    for (size_t i = 0; i < positions.size(); i++) {
        if (count[positions[i]] > 1) { preSkip[i] = 1; nrec++; }
    }
    if (nrec > 0) {
        LOG("Dropped " << nrec << " records at " << npos << " duplicated positions in " << source << " (all records at a duplicated position are removed).");
    }
}

// P1: CHROM, POS and ID of a kept VCF record, taken while the haplotypes
// are read, so the map does not need two more passes over the file. Fields
// are split on spaces/tabs as readMapDataVCF does; locId is the record's
// index among data lines, as there.
// M9: one reader for every VCF input. 3.1.0 had a single-population reader
// and a separate first/second pass for XP, each reading the file twice:
// once to check records and decide which to skip, once to load the alleles.
// readVCFRecords() reads the file once. Each record is checked, the
// per-file filters are applied (duplicated position, --filter-vcf, and MAF
// for single-population runs), and the alleles go into a scratch bitset
// that is kept only if the record is.
VCFRecords HapMap::readVCFRecords(const string &filename, bool xp, bool captureMap)
{
    VCFRecords rs;

    GzLineReader fin;   // P1: buffered reader
    fin.open(filename.c_str());
    if (fin.fail()) {
        HANDLE_ERROR("Failed to open " + filename + " for reading.");
    }

    constexpr int numMapCols = 9;   // fixed VCF columns before the genotype fields
    string line;

    int physpos = -1;
    int skip_due_to_duplicate_pos = 0;
    char lastSkipReason = 0;        // 'd', 'f' or 'm' for the last removed record
    bool prevKept = false;          // was the previous record kept?
    int nhaps = 0;
    MyBitset* scratchHap = nullptr;
    MyBitset* scratchXr = nullptr;

    auto find_field_end = [](const char* q, const char* end) -> const char* {
        while (q < end && *q != '\t') ++q;
        return q;
    };
    auto is_snp_base = [](const char* s, const char* e) -> bool {
        if (e - s != 1) return false;
        const char c = *s;
        return c == 'A' || c == 'C' || c == 'G' || c == 'T';
    };
    auto clearScratch = [&]() {
        std::memset(scratchHap->bits, 0, sizeof(uint64_t) * scratchHap->nwords);
        scratchHap->num_1s = 0;
        if (scratchXr) { std::memset(scratchXr->bits, 0, sizeof(uint64_t) * scratchXr->nwords); scratchXr->num_1s = 0; }
    };
    auto markSkipped = [&](size_t rec) { if (xp) rs.skipped[rec] = 1; rs.skipQueue.push(rec); };

    while (getline(fin, line))
    {
        if (line.empty()) continue;

        if (line[0] == '#') {
            if (line.rfind("#CHROM", 0) == 0) {
                int tabs = 0;
                for (char c : line) if (c == '\t') ++tabs;
                rs.current_ngts = (tabs + 1) - numMapCols;
                if (rs.current_ngts <= 0) {
                    HANDLE_ERROR("Failed to determine number of samples from VCF header in " + filename);
                }
                nhaps = p.UNPHASED ? rs.current_ngts : rs.current_ngts * 2;
                scratchHap = new MyBitset(nhaps);
                scratchXr = p.UNPHASED ? new MyBitset(nhaps) : nullptr;
            }
            continue;
        }
        if (rs.current_ngts < 0) {
            HANDLE_ERROR("VCF header (#CHROM line) not found before variant data in " + filename);
        }

        const size_t rec = rs.nloci_before_filtering++;
        if (xp) rs.skipped.push_back(0);

        const char* pcur = line.c_str();
        const char* end  = pcur + line.size();

        const char* field_starts[numMapCols];
        const char* field_ends[numMapCols];
        for (int col = 0; col < numMapCols; ++col) {
            field_starts[col] = pcur;
            field_ends[col]   = find_field_end(pcur, end);
            pcur = field_ends[col];
            if (pcur < end) ++pcur;
        }

        int new_physpos = 0;
        if (!parsePosition(field_starts[1], field_ends[1], new_physpos)) {   // B44
            HANDLE_ERROR("POS '" + string(field_starts[1], field_ends[1]) + "' of record " + to_string(rec + 1) + " in " + filename + POS_LIMIT_MSG);
        }

        if (xp) {
            // The XP reader pairs the two files by position, which needs both
            // sorted. Keep the position and a hash of REF<TAB>ALT per record.
            if (physpos >= 0 && new_physpos < physpos) {
                HANDLE_ERROR(filename + " is not sorted by position (" + to_string(new_physpos)
                             + " follows " + to_string(physpos) + ").");
            }
            rs.physpos.push_back((size_t)new_physpos);
            uint64_t h = 1469598103934665603ULL;            // FNV-1a 64
            for (const char* c = field_starts[3]; c < field_ends[4]; ++c) {
                h ^= (unsigned char)*c;
                h *= 1099511628211ULL;
            }
            rs.refalt_hash.push_back(h);
        }

        skip_due_to_duplicate_pos = (physpos == new_physpos) ? skip_due_to_duplicate_pos + 1 : 0;
        physpos = new_physpos;

        bool skip_due_to_missing = false;
        bool skip_due_to_multiallelic = false;
        if (p.MULTI_ALLELIC) {
            if (!is_snp_base(field_starts[3], field_ends[3]) || !is_snp_base(field_starts[4], field_ends[4])) {
                skip_due_to_multiallelic = true;
                rs.n_nonsnp++;
            }
        }

        // Genotypes (B16, B33): alleles other than 0/1 are an error (with
        // --filter-vcf the record is skipped instead), '/' is an error in
        // phased mode, each GT must be a diploid call, and the number of
        // genotype fields must match the #CHROM header.
        auto bad_site = [&](const string& msg) {
            if (!p.MULTI_ALLELIC) { HANDLE_ERROR(msg); }
            skip_due_to_missing = true;
        };

        int number_of_1s = 0, number_of_2s = 0;
        int seen_ngts = 0;
        uint64_t* const hapBits = scratchHap->bits;
        uint64_t* const xrBits = scratchXr ? scratchXr->bits : nullptr;
        while (pcur < end && seen_ngts < rs.current_ngts)
        {
            const char* g = pcur;
            const int shape = gtShape(g, end);
            if (shape == 0) {
                const char a1 = g[0];
                const char sep = g[1];
                const char a2 = g[2];
                if ((a1 != '0' && a1 != '1') || (a2 != '0' && a2 != '1')) {
                    bad_site("Alleles must be coded 0 or 1 only. Found alleles "
                             + string(1, a1) + sep + a2 + " at pos " + to_string(new_physpos)
                             + " in " + filename + ".");
                }
                if (!p.UNPHASED && sep != '|') {
                    bad_site("Unphased entries detected (/ is used) at pos " + to_string(new_physpos)
                             + " in " + filename + ". Make sure you run with --unphased flag for correct results.");
                }
                const int v1 = (a1 == '1');
                const int v2 = (a2 == '1');
                // bit indices are below nhaps by construction (seen_ngts <
                // current_ngts), so the bits are set without set_bit's check
                if (p.UNPHASED) {
                    const int sum = v1 + v2;
                    if (sum == 2) { number_of_2s++; xrBits[seen_ngts >> 6] |= 1ULL << (seen_ngts & 63); }
                    else if (sum == 1) { number_of_1s++; hapBits[seen_ngts >> 6] |= 1ULL << (seen_ngts & 63); }
                } else {
                    const int b = 2 * seen_ngts;
                    hapBits[b >> 6] |= ((uint64_t)v1 << (b & 63)) | ((uint64_t)v2 << ((b + 1) & 63));
                    number_of_1s += v1 + v2;
                }
                pcur = g + 3;   // the GT field is known up to here
            } else if (shape == 1) {
                bad_site("Haploid genotype '" + gtText(g, end) + "' at pos " + to_string(new_physpos) + " in " + filename
                         + ". selscan needs diploid genotypes (a|b, or a/b with " + ARG_UNPHASED + ").");
            } else {
                bad_site("Malformed genotype field '" + gtText(g, end) + "' at pos " + to_string(new_physpos) + " in " + filename + ".");
            }
            while (pcur < end && *pcur != '\t') ++pcur;
            if (pcur < end) ++pcur;
            seen_ngts++;
        }

        if (seen_ngts != rs.current_ngts) {
            HANDLE_ERROR("line " + to_string(rec + 1) + " of " + filename + " has " + to_string(seen_ngts) +
                         " genotype fields, but header indicates " + to_string(rs.current_ngts) + ".");
        }

        bool skip_due_to_maf = shouldSkipLocus(number_of_1s, number_of_2s, rs.current_ngts * 2);
        if (xp) skip_due_to_maf = false;   // XP is never MAF-filtered (Q1)

        // Every record at a duplicated position is removed. The first one is
        // removed when the second is seen: if it was kept, drop it now; if a
        // filter removed it, count it with the duplicates.
        if (skip_due_to_duplicate_pos == 1) {
            if (prevKept) {   // the previous record (rec - 1) was kept
                delete rs.hap.back(); rs.hap.pop_back();
                if (!rs.xr.empty()) { delete rs.xr.back(); rs.xr.pop_back(); }
                if (xp) rs.derived.pop_back();
                if (captureMap) rs.map.pop_back();
                markSkipped(rec - 1);
                rs.nSkipDup++;
            } else if (lastSkipReason == 'f' || lastSkipReason == 'm') {
                (lastSkipReason == 'f' ? rs.nSkipFilter : rs.nSkipMaf)--;
                rs.nSkipDup++;
            }
        }

        if (skip_due_to_maf || skip_due_to_multiallelic || skip_due_to_missing || skip_due_to_duplicate_pos) {
            markSkipped(rec);
            if (skip_due_to_duplicate_pos) { rs.nSkipDup++; lastSkipReason = 'd'; }
            else if (skip_due_to_multiallelic || skip_due_to_missing) { rs.nSkipFilter++; lastSkipReason = 'f'; }
            else { rs.nSkipMaf++; lastSkipReason = 'm'; }
            clearScratch();
            prevKept = false;
            continue;
        }

        // keep: the scratch bitsets become this locus's, and new ones are made
        scratchHap->num_1s = number_of_1s;
        rs.hap.push_back(scratchHap);
        if (scratchXr) { scratchXr->num_1s = number_of_2s; rs.xr.push_back(scratchXr); }
        prevKept = true;
        if (xp) rs.derived.push_back(number_of_1s + 2 * number_of_2s);   // unphased: hets + 2 * hom-derived
        if (captureMap) {
            MapEntry me;
            me.chr.assign(field_starts[0], field_ends[0]);
            me.physicalPos = new_physpos;
            me.locusName.assign(field_starts[2], field_ends[2]);
            me.locId = (int)rec;
            me.geneticPos = (long double)(new_physpos) / 1000000.0;
            rs.map.push_back(std::move(me));
        }
        scratchHap = new MyBitset(nhaps);
        scratchXr = p.UNPHASED ? new MyBitset(nhaps) : nullptr;
    }
    delete scratchHap;
    delete scratchXr;
    fin.close();
    return rs;
}

// Hand the kept loci to hapData, set its skip queue (record indices of the
// removed records, used by the --map reader) and, for the query file, the
// map entries taken from the VCF (P1).
void HapMap::installVCFRecords(VCFRecords &rs, HapData &hapData, bool captureMap)
{
    const int nhaps = p.UNPHASED ? rs.current_ngts : rs.current_ngts * 2;
    hapData.skipQueue = std::move(rs.skipQueue);
    vector<int>().swap(rs.derived);

    if (captureMap) {
        vcfMap.swap(rs.map);
        deque<MapEntry>().swap(rs.map);
    }
    hapData.adoptLoci(nhaps, rs.hap, rs.xr);   // takes the bitsets and frees both vectors
}


void HapMap::readHapDataVCFXP(string filename, string filename2, HapData& hapData, HapData& hapData2)
{
    initParamsInHap(hapData);
    initParamsInHap(hapData2);

    LOG("Reading " << filename << " and " << filename2 << ".");

    VCFRecords rs1 = readVCFRecords(filename, true, true);
    VCFRecords rs2 = readVCFRecords(filename2, true, false);
    // Pair the two files by position (B15). A record is kept only when its
    // position occurs exactly once in each file, neither record was removed
    // by the per-file filters (duplicate position, --filter-vcf), REF/ALT
    // agree, and the site is polymorphic in the two populations combined.
    // Kept records therefore appear in the same order in both files, which
    // the map reader (driven by file 1's skip queue) relies on.
    const size_t n1 = rs1.nloci_before_filtering;
    const size_t n2 = rs2.nloci_before_filtering;
    const vector<size_t>& pos1 = rs1.physpos;
    const vector<size_t>& pos2 = rs2.physpos;

    // Kept records are those not skipped, in file order, so k1 and k2 walk
    // the kept loci (bitsets and derived counts) alongside the records.
    vector<char> drop1(rs1.hap.size(), 0), drop2(rs2.hap.size(), 0);
    size_t only1 = 0, only2 = 0, partner1 = 0, partner2 = 0, mismatch = 0, mono = 0, shared = 0;
    const long total_alleles = 2L * (rs1.current_ngts + rs2.current_ngts);

    size_t i = 0, j = 0, k1 = 0, k2 = 0;
    while (i < n1 || j < n2) {
        if (j >= n2 || (i < n1 && pos1[i] < pos2[j])) {        // only in file 1
            if (!rs1.skipped[i]) { drop1[k1++] = 1; only1++; }
            i++;
            continue;
        }
        if (i >= n1 || pos2[j] < pos1[i]) {                     // only in file 2
            if (!rs2.skipped[j]) { drop2[k2++] = 1; only2++; }
            j++;
            continue;
        }
        // same position in both files
        const bool keep1 = !rs1.skipped[i], keep2 = !rs2.skipped[j];
        if (!keep1 || !keep2) {
            if (keep1) { drop1[k1] = 1; partner1++; }
            if (keep2) { drop2[k2] = 1; partner2++; }
        } else if (rs1.refalt_hash[i] != rs2.refalt_hash[j]) {
            drop1[k1] = 1; drop2[k2] = 1;
            mismatch++;
        } else {
            const long total_derived = (long)rs1.derived[k1] + rs2.derived[k2];
            if (total_derived == 0 || total_derived == total_alleles) {
                drop1[k1] = 1; drop2[k2] = 1;
                mono++;
            } else {
                shared++;
            }
        }
        if (keep1) k1++;
        if (keep2) k2++;
        i++;
        j++;
    }

    // remove the dropped loci, and rebuild each skip queue from the
    // per-file skips and the pairing drops
    auto applyDrops = [](VCFRecords &rs, const vector<char> &drop) {
        size_t w = 0, k = 0;
        const bool hasXr = !rs.xr.empty(), hasMap = !rs.map.empty();
        queue<size_t> q;
        for (size_t r = 0; r < rs.nloci_before_filtering; r++) {
            if (rs.skipped[r]) { q.push(r); continue; }
            if (drop[k]) {
                delete rs.hap[k]; if (hasXr) delete rs.xr[k];
                q.push(r);
            } else {
                if (w != k) {
                    rs.hap[w] = rs.hap[k]; rs.derived[w] = rs.derived[k];
                    if (hasXr) rs.xr[w] = rs.xr[k];
                    if (hasMap) rs.map[w] = std::move(rs.map[k]);
                }
                w++;
            }
            k++;
        }
        rs.hap.resize(w); rs.derived.resize(w);
        if (hasXr) rs.xr.resize(w);
        if (hasMap) rs.map.resize(w);
        rs.skipQueue = std::move(q);
    };
    applyDrops(rs1, drop1);
    applyDrops(rs2, drop2);
    vector<char>().swap(rs1.skipped); vector<char>().swap(rs2.skipped);
    vector<size_t>().swap(rs1.physpos); vector<size_t>().swap(rs2.physpos);
    vector<uint64_t>().swap(rs1.refalt_hash); vector<uint64_t>().swap(rs2.refalt_hash);

    if (only1 > 0)
        LOG("Dropped " << only1 << " records of " << filename << " at positions absent from " << filename2 << ".");
    if (only2 > 0)
        LOG("Dropped " << only2 << " records of " << filename2 << " at positions absent from " << filename << ".");
    if (partner1 > 0)
        LOG("Dropped " << partner1 << " records of " << filename << " whose matching record in " << filename2 << " was removed by a per-file filter (duplicated position, or --filter-vcf).");
    if (partner2 > 0)
        LOG("Dropped " << partner2 << " records of " << filename2 << " whose matching record in " << filename << " was removed by a per-file filter (duplicated position, or --filter-vcf).");
    if (mismatch > 0)
        LOG("Dropped " << mismatch << " shared positions where REF/ALT differ between the two files.");
    if (mono > 0)
        LOG("Dropped " << mono << " sites monomorphic in the two populations combined.");
    LOG(shared << " sites are shared by the two files and retained.");

    if (rs2.hap.empty()) {
        rs1.freeAll();
        HANDLE_ERROR("After filtering, no loci remain in reference VCF. Please check your input files."
                     << (rs2.n_nonsnp > 0 ? " " + ARG_MULTIALLELIC + " removed " + to_string(rs2.n_nonsnp)
                         + " records whose REF/ALT are not single A/C/G/T bases (for example 0/1 coding)." : string("")));
    } else if (rs1.hap.empty()) {
        rs2.freeAll();
        HANDLE_ERROR("After filtering, no loci remain in query VCF. Please check your input files."
                     << (rs1.n_nonsnp > 0 ? " " + ARG_MULTIALLELIC + " removed " + to_string(rs1.n_nonsnp)
                         + " records whose REF/ALT are not single A/C/G/T bases (for example 0/1 coding)." : string("")));
    }
    LOG("After filtering, " << rs1.hap.size() << " loci remain in query VCF.");
    LOG("After filtering, " << rs2.hap.size() << " loci remain in reference VCF.");
    LOG("XP mode: MAF-based filtering is disabled.");

    installVCFRecords(rs1, hapData, true);
    installVCFRecords(rs2, hapData2, false);
    LOG("====");
}

void HapMap::readHapDataVCF(string filename, HapData& hapData)
{
    initParamsInHap(hapData);
    LOG("Reading " << filename << ".");

    VCFRecords rs = readVCFRecords(filename, false, true);
    const size_t skipcount = rs.nloci_before_filtering - rs.hap.size();

    if (p.SKIP && !p.CALC_XP && !p.CALC_XPNSL) {
        LOG(ARG_SKIP << " set. Removing all variants < " << p.MAF << ".\n");
    }
    const int nhaps = p.UNPHASED ? rs.current_ngts : rs.current_ngts * 2;
    LOG("Loading " << nhaps << " haplotypes and " << rs.hap.size() << " loci. Skipped " << skipcount << " loci \n");

    if (rs.nloci_before_filtering == skipcount && rs.n_nonsnp > 0) {
        HANDLE_ERROR(ARG_MULTIALLELIC << " removed " << rs.n_nonsnp << " of " << rs.nloci_before_filtering
                     << " records of " << filename << " because REF/ALT are not single A/C/G/T bases"
                     << " (for example, files that code alleles as 0/1), and no loci remain.");
    }

    // 3.0.1 reported every removed record as a low frequency variant (M10)
    if (skipcount > 0) {
        string why;
        auto part = [&why](size_t n, const string &what) { if (n) why += (why.empty() ? "" : ", ") + to_string(n) + " " + what; };
        part(rs.nSkipDup, "at duplicated positions (all records at such a position)");
        part(rs.nSkipFilter, "removed by " + ARG_MULTIALLELIC);
        part(rs.nSkipMaf, "below the MAF cutoff");
        LOG("Removed " << skipcount << " records: " << why << ".");
    }

    installVCFRecords(rs, hapData, true);
    LOG("====");
}

// Decide whether to skip a locus based on MAF and missing data
bool HapMap::shouldSkipLocus(int number_of_1s, int number_of_2s, int nalleles_per_loc)
{
    // Compute the count of derived alleles depending on phased/unphased data
    int derived_allele_count = (p.UNPHASED
        ? (number_of_1s + number_of_2s * 2)  // unphased: 1s = het, 2s = homozygous derived
        : number_of_1s);                    // phased: count of 1s only

    // Flags for filtering
    bool skip_due_to_maf = false;     // Skip based on minor allele frequency threshold
    
    // Check for low-frequency alleles (MAF filtering)
    if (p.SKIP) {
        const long d = derived_allele_count, n = nalleles_per_loc;
        skip_due_to_maf = minor_allele_below(std::min(d, n - d), n, MIN_MAF_CUTOFF);
    }
    
    return skip_due_to_maf;

    /* UNCOMMENT if missing data is allowed
    bool skip_due_to_missing = false; // Skip based on missingness threshold
    if (p.MISSING_ALLOWED) {
        const double ALLOWED_MISSING_PERC = 1.0; // 100% allowed unless explicitly used
        int ancestral_allele_count = nalleles_per_loc - derived_allele_count;

        // Check if MAF threshold would be violated when accounting for missing data
        bool derived_too_rare = (derived_allele_count + missing_count) * 1.0 / nalleles_per_loc < MIN_MAF_CUTOFF;
        bool ancestral_too_rare = (ancestral_allele_count + missing_count) * 1.0 / nalleles_per_loc < MIN_MAF_CUTOFF;

        if (p.SKIP) {
            skip_due_to_missing = derived_too_rare || ancestral_too_rare;
        }

        // Check if missingness exceeds allowed threshold
        bool too_much_missing = (1.0 * missing_count / nalleles_per_loc) > ALLOWED_MISSING_PERC;
        skip_due_to_missing = skip_due_to_missing || too_much_missing;
    }
    return skip_due_to_maf || skip_due_to_missing;
    */


}


// bool HapMap::skipThisLocus(int number_of_1s, int number_of_2s, int nalleles_per_loc, int missing_count){ 
//     //nalleles same for phased and unphased, tot n of columns

//     int derived_allele_count = (p.UNPHASED? (number_of_1s + number_of_2s*2) : number_of_1s);
    
//     // if ((derived_allele_count)*1.0/(nalleles_per_loc) == 0){
//     //     return true;
//     // }

//     // if( 1-(derived_allele_count*1.0/(nalleles_per_loc)) == 0 ){
//     //     return true;
//     // }

//     bool skipreason1 = false; // --skip-low-freq filtering based on MAF
//     bool skipreason2 = false; // missing

//     double ALLOWED_MISSING_PERC = 1;
//     if(p.MISSING_ALLOWED){
//         int ancestral_allele_count = nalleles_per_loc - derived_allele_count;
//         bool derived_lower = (derived_allele_count+missing_count)*1.0/(nalleles_per_loc) < this->MIN_MAF_CUTOFF ;
//         bool ancestral_lower = (ancestral_allele_count+missing_count)*1.0/(nalleles_per_loc) < this->MIN_MAF_CUTOFF ;
//         if(p.SKIP){
//                 skipreason2 = (derived_lower || ancestral_lower);
//         }

//         bool missing_cutoff_crossed = 1.0*missing_count/nalleles_per_loc > ALLOWED_MISSING_PERC;
//         skipreason2 = (skipreason2 || missing_cutoff_crossed);
//         // if (missing_cutoff_crossed) {
//         //     //cout<<"Skipping cause Derived "<<derived_allele_count<<" Missing "<<missing_count<<" Total "<<current_nhaps*2<<endl;
//         //     skipreason2 = true;
//         // }
//     }
    
//     //if(!p.MISSING_ALLOWED){
//         skipreason1 = (p.SKIP && (derived_allele_count*1.0/(nalleles_per_loc) < this->MIN_MAF_CUTOFF || 1-(derived_allele_count*1.0/(nalleles_per_loc)) < this->MIN_MAF_CUTOFF ));    
//     //}
//     return (skipreason1 ||  skipreason2 );
// }





// TPED format description
// 1. Chromosome (string)
// 2. SNP ID (rs number: string)
// 3. Genetic distance (cM: float)
// 4. Physical position (int)
// 5. Allele 1: Sample 1 (0 or 1: char)
// 6. Allele 2: Sample 1 (0 or 1: char)
// 7. Allele 1: Sample 2 (0 or 1: char)
// 8. Allele 2: Sample 2 (0 or 1: char)
// ...
void HapMap::readHapDataTPED(string filename, HapData &hapData)
{
    initParamsInHap(hapData);
    igzstream fin;
    
    // PHASE1: opening file first time and counting number of loci and haps, and determine which loci to skip
    LOG("Opening " << filename << " to count number of haplotypes and loci...");

    fin.open(filename.c_str());
    if (fin.fail())
    {
        HANDLE_ERROR("Failed to open " + filename + " for reading.");
    }

    int numMapCols = 4;
    //int fileStart = fin.tellg();
    string line;
    int nloci_before_filter = 0; // number of loci before filtering (filtering happens when SKIP is enabled)
    int previous_nhaps = -1;
    int current_nhaps = 0;
    
    queue<int> skiplist;
    int n_maf_skipped = 0;

    //Counts number of haps (cols) and number of loci (rows)
    //if any lines differ, send an error message and throw an exception
    int lineNum = 0;
    while (getline(fin, line))
    {
        lineNum++;
        if (!hasData(line)) continue;   // U9
        nloci_before_filter++; // each line contributes to a locus
        current_nhaps = countFields(line) - numMapCols; // here, current_nhaps (count of fields) is the number of haps 
        
        stringstream ss;
        ss.str(line);

        string junk;
        
        for(int i = 0; i < numMapCols; i++){
            ss >> junk;
        }

        std::ostringstream rest;
        while (ss >> junk) {
            rest << junk ;
        }

        line.clear();
        string line_haps = rest.str();
        line_haps.erase(std::remove_if(line_haps.begin(), line_haps.end(), [](char c) {
            return std::isspace(static_cast<unsigned char>(c));
        }), line_haps.end());  // erase spaces inside and outside
        if (static_cast<int>(line_haps.size()) != current_nhaps) {
            HANDLE_ERROR("line " + std::to_string(lineNum) + " of " + filename + " has a genotype field longer than one character. Alleles must be coded 0/1 only.");
        }
        checkHapRow(line_haps, filename, lineNum, previous_nhaps < 0);

        int number_of_2s = 0; // for unphased
        int number_of_1s = 0;
        for (int hap = 0; hap < current_nhaps; hap++){ 
            if(p.UNPHASED){
                if (hap % 2 == 0){
                    if(line_haps[hap]=='1'  && line_haps[hap+1]=='1'){
                        number_of_2s++;
                    }else if ( (line_haps[hap]=='1' && line_haps[hap+1]=='0') || (line_haps[hap]=='0' && line_haps[hap+1]=='1') ){ // ==1
                        number_of_1s++;
                    }
                }
            }else{
                if(line_haps[hap]=='1'){
                    number_of_1s++;
                }
            }
        }
        if(preSkipAt(nloci_before_filter-1)){
            skiplist.push(nloci_before_filter-1);
        }else if(shouldSkipLocus(number_of_1s, number_of_2s, current_nhaps)){
            skiplist.push(nloci_before_filter-1);
            n_maf_skipped++;
        }

        if (previous_nhaps < 0)
        {
            previous_nhaps = current_nhaps;
            continue;
        }
        else if (previous_nhaps != current_nhaps) {
            HANDLE_ERROR(("line " + std::to_string(nloci_before_filter) +
                          " of " + filename + " has " + std::to_string(current_nhaps) +
                          " haps, but the previous line has " + std::to_string(previous_nhaps) + ".").c_str());
        }
        
        previous_nhaps = current_nhaps;
    }

    fin.clear(); // clear error flags
    //fin.seekg(fileStart);
    fin.close();


    if(p.SKIP) { //prefilter all sites < MAF
        LOG(ARG_SKIP << " set. Removing all variants < " << p.MAF << ".");
        LOG("Removed "<<n_maf_skipped <<" low frequency variants.");

    }else{
        LOG(ARG_KEEP << " set. NOT removing variants < " << p.MAF << ".");
    }
    
    if(p.UNPHASED){
        if (current_nhaps % 2 != 0)
        {
            HANDLE_ERROR("Number of haplotypes must be even for unphased.");
        }
    }
    current_nhaps = (p.UNPHASED) ? current_nhaps/2 : current_nhaps;
    
    LOG("Loading " << current_nhaps << " haplotypes and " << nloci_before_filter - skiplist.size() << " loci...");
    hapData.initHapData(current_nhaps, nloci_before_filter-skiplist.size());
    hapData.skipQueue = queue<size_t>();

    
    // PHASE 2: Open again and load the data to HapData: now that we gathered some idea about the format of the file
    fin.open(filename.c_str());

    if (fin.fail())
    {
        HANDLE_ERROR("Failed to open " + filename + " for reading.");
    }

    
    string junk;
    int locus_after_filter = 0;
    for (int locus = 0; locus < nloci_before_filter; locus++)
    {
        while (getline(fin, line) && !hasData(line)) {}   // U9
        if(!skiplist.empty()){
            if(skiplist.front() == locus){
                skiplist.pop();    
                hapData.skipQueue.push(locus); // make a copy to skipQueue
                continue;
            }
        }

        istringstream ss(line);
        for (int i = 0; i < numMapCols; i++)
        {
            ss >> junk; // skip first 4 columns
        }

        std::ostringstream rest;

        int cols = current_nhaps; 
        
        if(p.UNPHASED) {
            cols *= 2;
        }
        for (int i = 0; i < cols; i++)
        {
            ss >> junk; 
            rest << junk << "";
        }

        string line_haps = rest.str();
        
        line_haps.erase(std::remove_if(line_haps.begin(), line_haps.end(), [](char c) {
            return std::isspace(static_cast<unsigned char>(c));
        }), line_haps.end());  // erase spaces inside and outside

        if (p.UNPHASED){
            for (int hap = 0; hap < cols; hap++){ 
                if (hap % 2 == 0){
                    if(line_haps[hap]=='1'  && line_haps[hap+1]=='1'){
                        hapData.addAllele2(locus_after_filter, hap/2);
                    }else if ( (line_haps[hap]=='1' && line_haps[hap+1]=='0') || (line_haps[hap]=='0' && line_haps[hap+1]=='1') ){ // ==1
                        hapData.addAllele1(locus_after_filter, hap/2);
                    }
                }
            }
        }else{ // PHASED (0/1 coding checked in pass 1)
            for (int hap = 0; hap < current_nhaps; hap++)
            {
                char allele1 = line_haps[hap];
                if(allele1=='1'){
                    hapData.addAllele1(locus_after_filter, hap);
                }
            }
        }
        locus_after_filter++;
        line.clear();
        ss.clear();
        rest.clear();
    }

    fin.close();
}

void HapMap::loadHapMapData(){
        
    mapData = std::make_unique<MapData>(); 
    hapData = std::make_unique<HapData>();

    // PHASE 2:  Load data to HapData and MapData
    if (p.CALC_XP || p.CALC_XPNSL){
        hapData2 = std::make_unique<HapData>();
    }   

    auto start_reading = std::chrono::high_resolution_clock::now();

    // Sites at duplicated positions are dropped by every reader (preSkip).
    preSkip.clear();
    if (p.TPED) {
        vector<long> pos = readPositionColumn(p.tpedFilename);
        if (p.CALC_XP || p.CALC_XPNSL) {
            if (readPositionColumn(p.tpedFilename2) != pos) {
                HANDLE_ERROR(p.tpedFilename + " and " + p.tpedFilename2 + " do not list the same positions in the same order.");
            }
        }
        markDuplicatePositions(pos, p.tpedFilename);
    } else if (!p.VCF) {
        markDuplicatePositions(readPositionColumn(p.mapFilename), p.mapFilename);
    }

    // XP: drop sites monomorphic in the two populations combined, as the VCF
    // XP reader does. TPED and hap/thap kept them, and the extra sites change
    // iHH through the integration steps and gap scaling, so the same data gave
    // different XP results by format.
    if ((p.CALC_XP || p.CALC_XPNSL) && !p.VCF) {
        const char fmt = p.TPED ? 't' : (p.THAP ? 'h' : 'm');
        const string f1 = p.TPED ? p.tpedFilename : (p.THAP ? p.thapFilename : p.hapFilename);
        const string f2 = p.TPED ? p.tpedFilename2 : (p.THAP ? p.thapFilename2 : p.hapFilename2);
        vector<long> d1, d2; long n1 = 0, n2 = 0;
        scanDerivedCounts(f1, fmt, d1, n1);
        scanDerivedCounts(f2, fmt, d2, n2);
        if (d1.size() != d2.size()) {
            HANDLE_ERROR(f1 + " and " + f2 + " do not have the same number of loci.");
        }
        if (preSkip.size() < d1.size()) preSkip.resize(d1.size(), 0);
        size_t nmono = 0;
        for (size_t i = 0; i < d1.size(); i++) {
            const long d = d1[i] + d2[i];
            if (!preSkip[i] && (d == 0 || d == n1 + n2)) { preSkip[i] = 1; nmono++; }
        }
        if (nmono > 0) {
            LOG("Dropped " << nmono << " sites monomorphic in the two populations combined.");
        }
    }

    if (p.TPED)
    {
        readHapDataTPED(p.tpedFilename, *hapData);
        if (p.CALC_XP || p.CALC_XPNSL)
        {
            readHapDataTPED(p.tpedFilename2, *hapData2);
            if (hapData->nloci != hapData2->nloci)
            {
                HANDLE_ERROR("Haplotypes from " + p.tpedFilename + " and " + p.tpedFilename2 + " do not have the same number of loci.");
            }
        }
        int exp_haps = (p.UNPHASED)? hapData->nhaps*2 : hapData->nhaps;
        mapData->readMapDataTPED(p.tpedFilename, hapData->nloci, exp_haps, p.USE_PMAP, hapData->skipQueue);
    }
    else if (p.VCF) {
        if (p.CALC_XP || p.CALC_XPNSL)
        {
            // readHapDataVCF(p.vcfFilename, *hapData); //query_vcf
            // readHapDataVCF(p.vcfFilename2, *hapData2); //ref_vcf
            readHapDataVCFXP(p.vcfFilename, p.vcfFilename2,  *hapData, *hapData2); // read both vcfs together and determine skiplists based on both
            // if (hapData->nloci != hapData2->nloci)
            // {
            //     HANDLE_ERROR("Haplotypes from " + p.vcfFilename + " and " + p.vcfFilename2 + " do not have the same number of loci.");
            // }
        }else{
            readHapDataVCF(p.vcfFilename, *hapData);
        }
        // B50: use --map whenever one is given (and --pmap is not set),
        // whatever statistics are requested. 3.0.1 ignored the map when
        // --nsl or --xpnsl was also requested, so iHS/XP-EHH in such a
        // combined run used physical distance.
        if(!p.USE_PMAP && p.mapFilename.compare(DEFAULT_FILENAME_MAP) != 0) {
            mapData->readMapData(p.mapFilename, hapData->nloci, p.USE_PMAP, hapData->skipQueue);
            // U8: map rows are matched to VCF records by order, so they must
            // agree on position (the output takes chr, id and pos from the
            // map). 3.0.1 compared only the counts.
            if (vcfMap.size() == (size_t)mapData->nloci) {
                for (int i = 0; i < mapData->nloci; i++) {
                    const MapEntry &m = mapData->mapEntries[i];
                    if (m.physicalPos != vcfMap[i].physicalPos) {
                        HANDLE_ERROR("data row " + to_string(m.locId + 1) + " of " + p.mapFilename + " (" + m.locusName + ") has position "
                                     + to_string(m.physicalPos) + ", but the matching record of " + p.vcfFilename + " has POS "
                                     + to_string(vcfMap[i].physicalPos) + ". With " + ARG_FILENAME_MAP + ", the map must list the VCF records in the same order, at the same positions.");
                    }
                }
            }
        }
        else if (vcfMap.size() == (size_t)hapData->nloci) {
            // P1: physical positions recorded while the haplotypes were read
            mapData->setFromEntries(vcfMap);
        }
        else{//Load physical positions
            mapData->readMapDataVCF(p.vcfFilename, hapData->nloci, hapData->skipQueue);
        }
        deque<MapEntry>().swap(vcfMap);
    }else if(p.THAP){
        readHapDataTHAP(p.thapFilename, *hapData);
        if (p.CALC_XP || p.CALC_XPNSL)
        {
            readHapDataTHAP(p.thapFilename2, *hapData2);
            if (hapData->nloci != hapData2->nloci)
            {
                HANDLE_ERROR("Haplotypes from " + p.thapFilename + " and " + p.thapFilename2 + " do not have the same number of loci.");
            }
        }
        mapData->readMapData(p.mapFilename, hapData->nloci, p.USE_PMAP, hapData->skipQueue);
    }else
    {
        readHapDataMSHap(p.hapFilename, *hapData);
        if (p.CALC_XP || p.CALC_XPNSL)
        {
            readHapDataMSHap(p.hapFilename2, *hapData2);
            if (hapData->nloci != hapData2->nloci)
            {
                HANDLE_ERROR("Haplotypes from " + p.hapFilename + " and " + p.hapFilename2 + " do not have the same number of loci.");
            }
        }
        mapData->readMapData(p.mapFilename, hapData->nloci, p.USE_PMAP, hapData->skipQueue);
    }

    // PHASE 2:  XOR
    // if(p.benchmark_flag=="XOR")
    //     hapData->xor_for_phased_and_unphased();

    auto end_reading = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> read_duration =  end_reading - start_reading;
    LOG("Input file loaded in "+to_string(read_duration.count())+" s.");
    LOG("=====");

    // DEBUG::: mapData.print();

    // if(mapData->chr_list.size() > 1 && p.MULTI_CHR == false){
    //     cerr<<"ERROR: Input file contains multiple chromosomes, but no chromosome list is provided. Use --multi-chr flag to specify chromosomes.\n";
    // }

    // Check if map is in order
    for (int i = 1; i < mapData->nloci; i++) {
        if ( mapData->mapEntries[i].physicalPos <  mapData->mapEntries[i-1].physicalPos ) {
            string ERROR_MSG = "Variant physical position must be monotonically increasing.\n";
            ERROR_MSG += "\t" + mapData->mapEntries[i].locusName + " " + std::to_string(mapData->mapEntries[i].physicalPos) + " appears after";
            ERROR_MSG += "\t" + mapData->mapEntries[i-1].locusName + " " + std::to_string(mapData->mapEntries[i-1].physicalPos) + ".";
            HANDLE_ERROR(ERROR_MSG);
        }
        if ( !p.CALC_NSL && mapData->mapEntries[i].geneticPos  < mapData->mapEntries[i-1].geneticPos  ) {
            string ERROR_MSG = "Variant genetic position must be monotonically increasing.\n";
            ERROR_MSG += "\t" + mapData->mapEntries[i].locusName + " " + std::to_string(mapData->mapEntries[i].geneticPos) + " appears after";
            ERROR_MSG += "\t" + mapData->mapEntries[i-1].locusName + " " + std::to_string(mapData->mapEntries[i-1].geneticPos) + ".";
            HANDLE_ERROR(ERROR_MSG);
        }
    }
    
}