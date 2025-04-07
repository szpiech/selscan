# selscan -- a program to calculate EHH-based scans for positive selection in genomes

Copyright (C) 2014  Zachary A Szpiech

selscan currently implements EHH, iHS, XP-EHH, nSL, XP-nSL and iHH12.    

It should be run separately for each chromosome and population (or population    
pair for XP-EHH).  selscan is 'dumb' with respect ancestral/derived coding and   
simply expects haplotype data to be coded 0/1.  Unstandardized iHS/nSL scores are    
thus reported as log(iHH1/iHH0) based on the coding you have provided.   


## 🛠️ Installation from source 

```bash
git clone https://github.com/szpiech/selscan/
cd selscan && git checkout speedup
cd src && make -f Makefile_platform

Replace Makefile_platform with:

Makefile for macOS
Makefile_linux for Linux
Makefile_win for Windows

### 📦 Precompiled Binaries

Precompiled binaries are available for the following platforms:

- **macOS (Apple Silicon, ARM64):** `/bin/macos/`  
- **macOS (Intel, x86_64):** `/bin/osx/`  
- **Windows:** `/bin/win/`  
- **Linux:** `/bin/linux/`


## Usage   
```
** Data must have no missing genotypes. **

selscan v2.0.0 -- a program to calculate EHH-based scans for positive selection in genomes.
Source code and binaries can be found at <https://www.github.com/szpiech/selscan>.

selscan currently implements EHH, iHS, XP-EHH, nSL, and XP-nSL.

To calculate EHH:

./selscan --ehh <locusID> --vcf <vcf> --map <mapfile> --out <outfile>

To calculate iHS:

./selscan --ihs --vcf <vcf> --map <mapfile> --out <outfile>

To calculate nSL:

./selscan --nsl --vcf <vcf> --out <outfile>

To calculate XP-nSL:

./selscan --xpnsl --vcf <vcf> --vcf-ref <vcf> --out <outfile>

To calculate iHH12:

./selscan --ihh12 --vcf <vcf> --map <mapfile> --out <outfile>

To calculate XP-EHH:

./selscan --xpehh --vcf <vcf> --vcf-ref <vcf> --map <mapfile> --out <outfile>

----------Command Line Arguments----------

--alt <bool>: Set this flag to calculate homozygosity based on the sum of the
	squared haplotype frequencies in the observed data instead of using
	binomial coefficients.
	Default: false

--cutoff <double>: The EHH decay cutoff.
	Default: 0.05

--ehh <string>: Calculate EHH of the '1' and '0' haplotypes at the specified
	locus. Output: <physical dist> <genetic dist> <'1' EHH> <'0' EHH>
	Default: __NO_LOCUS__

--ehh-win <int>: When calculating EHH, this is the length of the window in bp
	in each direction from the query locus.
	Default: 100000

--gap-scale <int>: Gap scale parameter in bp. If a gap is encountered between
	two snps > GAP_SCALE and < MAX_GAP, then the genetic distance is
	scaled by GAP_SCALE/GAP.
	Default: 20000

--hap <string>: A hapfile with one column per haplotype, and one row per
	variant. Variants should be coded 0/1
	Default: __hapfile1

--help <bool>: Prints this help dialog.
	Default: false

--ihh12 <bool>: Set this flag to calculate iHH12.
	Default: false

--ihs <bool>: Set this flag to calculate iHS.
	Default: false

--ihs-detail <bool>: Set this flag to write out left and right iHH scores for '1' and '0' in addition to iHS.
	Default: false

--keep-low-freq <bool>: Include low frequency variants in the construction of your haplotypes.
	Default: false

--maf <double>: If a site has a MAF below this value, the program will not use
	it as a core snp.
	Default: 0.05

--map <string>: A mapfile with one row per variant site.
	Formatted <chr#> <locusID> <genetic pos> <physical pos>.
	Default: __mapfile

--max-extend <int>: The maximum distance an EHH decay curve is allowed to extend from the core.
	Set <= 0 for no restriction.
	Default: 1000000

--max-extend-nsl <int>: The maximum distance an nSL haplotype is allowed to extend from the core.
	Set <= 0 for no restriction.
	Default: 100

--max-gap <int>: Maximum allowed gap in bp between two snps.
	Default: 200000

--nsl <bool>: Set this flag to calculate nSL.
	Default: false

--out <string>: The basename for all output files.
	Default: outfile

--pi <bool>: Set this flag to calculate mean pairwise sequence difference in a sliding window.
	Default: false

--pi-win <int>: Sliding window size in bp for calculating pi.
	Default: 100

--pmap <bool>: Use physical map instead of a genetic map.
	Default: false

--ref <string>: A hapfile with one row per haplotype, and one column per
	variant. Variants should be coded 0/1. This is the 'reference'
	population for XP-EHH calculations.  Ignored otherwise.
	Default: __hapfile2

--skip-low-freq <bool>: **This flag is now on by default. If you want to include low frequency variants
in the construction of your haplotypes please use the --keep-low-freq flag.
	Default: false

--threads <int>: The number of threads to spawn during the calculation.
	Partitions loci across threads.
	Default: 1

--tped <string>: A TPED file containing haplotype and map data.
	Variants should be coded 0/1
	Default: __hapfile1

--tped-ref <string>: A TPED file containing haplotype and map data.
	Variants should be coded 0/1. This is the 'reference'
	population for XP-EHH calculations and should contain the same number
	of loci as the query population. Ignored otherwise.
	Default: __hapfile2

--trunc-ok <bool>: If an EHH decay reaches the end of a sequence before reaching the cutoff,
	integrate the curve anyway (iHS and XPEHH only).
	Normal function is to disregard the score for that core.
	Default: false

--unphased <bool>: Set this flag to use multilocus genotypes.
	Default: false

--vcf <string>: A VCF file containing haplotype data.
	A map file must be specified with --map.
	Default: __hapfile1

--vcf-ref <string>: A VCF file containing haplotype and map data.
	Variants should be coded 0/1. This is the 'reference'
	population for XP-EHH calculations and should contain the same number
	of loci as the query population. Ignored otherwise.
	Default: __hapfile2

--wagh <bool>: Set this flag to calculate XP-EHH using definition of EHH which
	separates core SNP alleles in the denominator.
	Default: false

--xpehh <bool>: Set this flag to calculate XP-EHH.
	Default: false

--xpnsl <bool>: Set this flag to calculate XP-nSL.
	Default: false
```

## Citations
```
ZA Szpiech (2021) selscan 2.0: scanning for sweeps in unphased data. biorxiv doi: 
	doi:10.1101/2021.10.22.465497.
ZA Szpiech and RD Hernandez (2014) selscan: an efficient multi-threaded program 
	to calculate EHH-based scans for positive selection. Molecular Biology and Evolution 
	31: 2824-2827.
ZA Szpiech et al. (2021) Application of a novel haplotype-based scan for local adaptation 
	to study high-altitude adaptation in rhesus macaques. Evolution Letters 
	doi: https://doi.org/10.1002/evl3.232
R Torres et al. (2018) Human demographic history has amplified the effects of
	background selection across the genome. PLoS Genetics 15: e1007898.
N Garud et al. (2015) Recent selective sweeps in North American Drosophila
	melanogaster show signatures of soft sweeps. PLoS Genetics 11: 1–32.
A Ferrer-Admetlla et al. (2014) On detecting incomplete soft or hard selective sweeps
	using haplotype structure. Molecular Biology and Evolution 31: 1275-1291.
K Wagh et al. (2012) Lactase Persistence and Lipid Pathway Selection in the Maasai. PloS ONE 7: e44751.
PC Sabeti et al. (2007) Genome-wide detection and characterization of positive 
	selection in human populations. Nature 449: 913–918.
BF Voight et al. (2006) A map of recent positive selection in the human 
	genome. PLoS Biology 4: e72.
PC Sabeti et al. (2002) Detecting recent positive selection in the human 
	genome from haplotype structure. Nature 419: 832–837.
```


## Change Log
```
17NOV2023 - selscan v2.0.1 - Bug fixes for --ehh flag. Total EHH at the core snp will now be reported correctly (i.e. homozygosity of the site and not as 0). Also implemented --unphased for --ehh, and EHH output files now have a header line.

	   - selscan v2.0.2 - Small change that should result in faster runtime when --pmap set.

22OCT2021 - selscan v2.0.0 - Introducing unphased versions of iHS, nSL, XP-EHH, and XP-nSL. Use with --unphased flag. See ZA Szpiech (2021) Biorxiv for details. Normalize as you would with the phased 
	statistics.
	
20MAY2020 - selscan v1.3.0 - Log ratios are now output as log10 not natural logs (beware comparisons with raw selscan computations from versions prior to v1.3.0). New statistics implemented.

	--pmap <bool>: Set this flag to use physical distance instead of genetic map

Introduction of XP-nSL, this statistic is a cross population statistic for identifying hard/soft sweeps. Does not require a genetic map. XP-nSL:nSL::XP-EHH:iHS

	--xpnsl <bool>: Set this flag to calculate XP-nSL.
	Default: false

Normalize XP-nSL with --xpnsl flag in norm.

lasugden adds the option to calculate XP-EHH with either definition of EHH. By default, uses original denominator (N choose 2). To use denominator defined in Wagh et al. for better performance on incomplete sweeps, use flag --wagh

	--wagh <bool>: Set this flag to calculate EHH with Wagh denominator. For xpehh only. DO NOT use with --alt
	Default: false

Normalize these computations with --xpehh flag in norm.

norm v1.3.0 - Now supports --xpnsl flag, which is identical to using --xpehh.
--qbins now has a default value of 10 instead of 20.
--bp-win analyses have been changed when analyzing XP-EHH and XP-nSL scores. Since positive scores suggest adaptation in the first (non-ref) population and negative scores suggest adaptation in the second (ref) population, we split windows into those enriched for extreme positive scores and those enriched for extreme negative scores.
min and max scores are given for each window for XP statistics, and the max |score| is reported for iHS and nSL stats.

*.windows output files therefore have additional columns:

For XP stats:
<win start> <win end> <# scores in win> <frac scores gt threshold> <frac scores lt threshold> <approx percentile for gt threshold wins> <approx percentile for lt threshold wins> <max score> <min score>

For iHS and nSL:
<win start> <win end> <# scores in win> <frac scores gt threshold> <frac scores lt threshold> <approx percentile for gt threshold wins> <approx percentile for lt threshold wins> <max score> <min score>


18SEPT2017 - norm v1.2.1a, selcan v1.2.0a iHH12 output files have a header line, XPEHH header line has an extra column name, fixed norm bugs relating to normalization of ihh12 files.

25AUG2017 - norm v1.2.1 released to fix a crash when --nsl flag is used.

18JUL2017 - Support for iHH12 calculations. norm has --ihh12 and --nsl flags.

09JAN2017 - Fixed buggy --crit-percent flag in norm binary.

05SEPT2016 - Fixed misleading error messages when --trunc-ok used.

12FEB2016 - v 1.1.0b - The flag --skip-low-freq is now on by default and no longer has any function.  selscan now filters low frequency variants by default.  A new flag --keep-low-freq is available if you would like to include low frequency variants when building haplotypes (low frequency variants will still be skipped over as core loci), using this option may reduce the power of iHS scans.

28OCT2015 - Updates to norm so that it can handle output from selscan when --ihs-detail is used.

18JUNE2015 - v1.1.0a - When calculating nSL, a mapfile is no longer required for VCF.  Physical distances will be read directly from VCF.  A mapfile specifying physical distances is stille required for .hap files when calculating nSL.  selscan now appropriately reports an error if this is not provided.

15JUNE2015 - Release of 1.1.0. tomkinsc adds the --ihs-detail parameter which, when provided as an adjunct to --ihs, will cause selscan to write out four additional columns to the output file of iHS calculations (in order): derived_ihh_left, derived_ihh_right, ancestral_ihh_left, and ancestral_ihh_right.

An example file row follows, with header added for clarity.

locus           phys-pos        1_freq          ihh_1           ihh_0           ihs             derived_ihh_left     derived_ihh_right    ancestral_ihh_left      ancestral_ihh_right

16133705        16133705        0.873626        0.0961264       0.105545        -0.0934761      0.0505176             0.0456087           0.0539295               0.0516158

From these values we can calculate iHS, but it is preserved in the output for convenience. Having left and right integral information may assist certain machine learning models that gain information from iHH asymmetry. 

selscan can now calculate the nSL statistic described in A Ferrer-Admetlla, et al. (2014) MBE, 31: 1275-1291.  Also introduced a check on map distance ordering.  Three new command line options.

--nsl <bool>: Set this flag to calculate nSL.
	Default: false

--max-extend-nsl <int>: The maximum distance an nSL haplotype is allowed to extend from the core.
	Set <= 0 for no restriction.
	Default: 100

--ihs-detail <bool> : Set this flag to write out left and right iHH scores for '1' and '0' in addition to iHS.

06MAY2015 - Release of 1.0.5. Added basic VCF support.  selscan can now read .vcf and .vcf.gz files but without tabix support.  A mapfile is required when using VCF.  Two new command line options.

13MAY2015 - norm v1.0.5 is released.  norm will now normalize ihs or xpehh scores.  Two new command line options.

--ihs <bool>: Do iHS normalization.

--xpehh <bool>: Do XP-EHH normalization.

Exactly one of these must be specified when running norm (e.g. ./norm --ihs --files *.ihs.out or ./norm --xpehh --files *.xpehh.out).

--vcf <string>: A VCF file containing haplotype data.
	A map file must be specified with --map.

--vcf-ref <string>: A VCF file containing haplotype and map data.
	Variants should be coded 0/1. This is the 'reference'
	population for XP-EHH calculations and should contain the same number
	of loci as the query population. Ignored otherwise.

07JAN2015 - norm bug fix and --skip-low-freq works for single EHH queries.

12NOV2014 - The program norm has been updated to allow for user defined critical values.  Two new command line options.

--crit-percent <double>: Set the critical value such that a SNP with iHS in the most extreme CRIT_PERCENT tails (two-tailed) is marked as an extreme SNP.
	Not used by default.

--crit-val <double>: Set the critical value such that a SNP with |iHS| > CRIT_VAL is marked as an extreme SNP.  Default as in Voight et al.
	Default: 2.00

17OCT2014 - Release of 1.0.4. A pairwise sequence difference module has been introduced.  This module isn't multithreaded at the moment, but still runs quite fast.  Calculating pi in 100bp windows with 198 haplotypes with 707,980 variants on human chr22 finishes in 77s on the test machine.  Using 100kb windows, it finishes in 34s.  Two new command line options.

--pi <bool>: Set this flag to calculate mean pairwise sequence difference in a sliding window.
	Default: false

--pi-win <int>: Sliding window size in bp for calculating pi.
	Default: 100

15SEP2014 - Release of 1.0.3.  **A critical bug in the XP-EHH module was introduced in version 1.0.2 and had been fixed in 1.0.3.  Do not use 1.0.2 for calculating XP-EHH scores.**  Thanks to David McWilliams for finding this error.  1.0.3 also introduces support for gzipped input files.  You may pass hap.gz, map.gz. and tped.gz files interchangably with unzipped files using the same command line arguments.  A new command line option is available.

--trunc-ok <bool>: If an EHH decay reaches the end of a sequence before reaching the cutoff,
	integrate the curve anyway (iHS and XPEHH only).
	Normal function is to disregard the score for that core.
	Default: false

17JUN2014 - Release of 1.0.2.  General speed improvements have been made, especially with threading.  New support for TPED formatted data and new command line options are available.

--skip-low-freq <bool>: Do not include low frequency variants in the construction of haplotypes (iHS only).
	Default: false

--max-extend: The maximum distance an EHH decay curve is allowed to extend from the core.
	Set <= 0 for no restriction.
	Default: 1000000

--tped <string>: A TPED file containing haplotype and map data.
	Variants should be coded 0/1
	Default: __hapfile1

--tped-ref <string>: A TPED file containing haplotype and map data.
	Variants should be coded 0/1. This is the 'reference'
	population for XP-EHH calculations and should contain the same number
	of loci as the query population. Ignored otherwise.
	Default: __hapfile2

10APR2014 - Release of 1.0.1.  Minor bug fixes. XP-EHH output header is now separated by tabs instead of spaces.  Removed references to missing data (which is not accepted), and introduced error checking in the event of non-0/1 data being provided.

26MAR2014 - Initial release of selscan 1.0.0.
```
