# Branch `code-review`: review fixes on top of 3.0.1

Base: `main` at 7a93ff1 ("Patch 3_0_1"). Finding ids (B1, U12, ...) refer to
the code review report (selscan_code_review.md). Every commit passes
`make check`. Since 341bdc0 (M1) the macOS and Linux Makefiles track header
dependencies, so an incremental build is safe. Before that commit, run
`make clean` after editing a header.

## Commits and what they change in the output

| commit | finding | output change |
|---|---|---|
| ef480b7 | harness | adds tests/ and `make check`; 30 golden tables recorded from 7a93ff1 |
| 78b076b | B5 | none (`--ehh`/`--ehh12` no longer exit 139 after writing output) |
| b47d8c2 | U12 | none (a missing `--ehh12` locus exits 1) |
| 7bf4094 | B3 | TPED XP now reads `--tped-ref`; byte-identical to the same data as VCF |
| a83316a | B1 | phased iHS/nSL drop cores that must cross a gap > `--max-gap`; gap fixture only (768 iHS, 200 nSL cores); all remaining rows identical |
| 8e701b8 | B2 | XP (VCF) restores sites at joint frequency 0.5 (14 phased) and sites without heterozygotes (9 unphased); nearby values move |
| feec363 | B6 | unphased XP-EHH/XP-nSL population-2 columns change (493/770 rows); `--wagh` p2 normaliser typo |
| d45021f | B4 | norm: XP/iHH12 p1 = 1.0 sites normalised in the single bin; iHS/nSL freq = 1.0 excluded and reported; norm-xpehh only |
| 3205f67 | B10 | unphased XP no longer skips sites lacking one homozygous class (770 -> 2193 sites; existing rows identical) |
| cd3c6f2 | B11 | unphased iHS cutoff with homozygous classes of < 2 individuals (1416 -> 1419 cores; 254 ihs values change) |
| 2903bfe | notes | this file |
| 341bdc0 | M1 | none (build: header dependencies via -MMD -MP) |
| f97de30 | B12 | none in goldens; `--pi` correct at >= 46,342 haplotypes (was negative) |
| 3558905 | B13 | none in goldens; `--pi` no longer reads past the last site when it is 0 or -1 mod `--pi-win` |
| fcea18f | B7 | phased nSL integrates with no EHH cutoff; nsl, nsl-thap, nsl-gap, norm-nsl change |
| 7b54c0e | tests | adds a 100 kb gap fixture (ihs/nsl/ihh12-gap100), the first golden to exercise `--gap-scale` |
| 02a4289 | B8 | iHH12 applies the gap scale once; ihh12-gap100 only (266 of 5477 cores) |
| 1378603 | Q5 | iHH12 truncates on EHH12 (definition change); every iHH12 value rises, 1361 -> 1145 cores on p1 (more reach a chromosome end) |
| 4de29fa | B9 | `--alt` uses squared frequencies (helper for iHS/nSL/EHH/unphased); ihs-alt changes in all rows |
| 28de5ec | B9 | `--ehh12 --alt` output normalised by N^2; EHH2/1 under `--alt`; new ehh12-alt golden |
| 13ffd40 | new | EHH2/1 at the core removes the larger class (was the derived class); one cell of ehh12; new ehh12-minor golden |
| b69fead | B41 | none (square_alt in double) |
| 0fc9333 | B27 | none in goldens; `--keep-low-freq` honoured in combined and `--multi-param` runs |
| 82d0042 | B32 | none in goldens; one symmetric MAF test with ties kept; at `--maf 0.1` on p1 iHS goes 178 -> 198 cores |
| a1c5bbd | B17 | none in goldens; `--pi`/`--ehh` combined with iHS see all sites |
| f985491 | data | none in goldens; TPED and hap/thap drop every record at a duplicated position, as VCF did, and log it |
| 28f0f22 | data | none in goldens; TPED/thap/hap XP drop sites monomorphic in the two populations combined, as VCF did (1663 on p1/p2) |
| 56de517 | B15 | none in goldens; XP VCF files are paired by position and REF/ALT (was by record index); unpaired and mismatched records are dropped and counted |
| 22769f4 | B47 | XP-nSL window and unphased step counted in retained sites (were raw file-1 records); xpnsl and xpnsl-unphased change in every row |
| aba04c4 | B16 | none in goldens; XP VCF reader rejects missing alleles, '/' in phased mode, short rows and a missing #CHROM line, as the single-population reader does |
| 5bd30f8 | tests | comment states measured 3.0.1 behaviour for the xp-validate fixtures |
| a018ca0 | B23 | none in goldens; `--norm-files --gene-bed` no longer crashes after writing |
| d2bd8ed | B24 | none in goldens; `--first` with several `--files` and `--gene-bed` stays in bounds |
| 9fbb6d7 | B25 | none in goldens; gene table with fewer than 3 genes: NA in `*_adj`, no GSL abort |
| 4535177 | B26 | none in goldens; XP `bottom_score_adj` is the residual of the bottom score |
| 49d0827 | B37 | none in goldens; a gene's last base counts for SNP and window overlap |
| 9f30ab1 | B48 | none in goldens; `--bp-win --gene-bed` annotates the `.windows` table (was the `.norm` file) |
| 1156907 | B38 | none in goldens; permutation test reads `top_score_adj`, add-one p-value, `--perm-seed` |
| 3078566 | B18 | none in goldens; non-numeric/non-finite scores are treated as missing and counted |
| 9d50efc | B34 | frequency bins exactly left-closed; norm-ihs, norm-nsl, norm-ihs-bpwin (+windows) re-recorded |
| a53b5b2 | B19, U10 | none in goldens; `--crit-percent` from finite normalised scores; iHH12 one-tailed |
| 4545108 | B21 | windows golden gains its last window (5 SNPs, below `--min-snps`) |
| ff34218 | B45 | none in goldens; iHH12 window maxima are per file |
| 3213f74 | B22, U4 | windows golden: 8 labels 1 -> 100; quantile bins sliced with the output rule, `>` rule, warning below 100 windows per bin |
| c3c2060 | B49 | none in goldens; norm `--log` written in every mode and kept open to the end |
| 15426a0 | B20 | none in goldens; `--log-input` takes the bin count from the log; bin statistics logged at 17 digits |
| 8150dd1 | U2 | none in goldens; `--log-only` honoured for iHS/nSL, window and gene analyses skipped |
| 02aaea8 | B39 | none in goldens; errors exit 1 with a message (were exit 134); `--multi-param` keys and types checked |
| 1ff78da | U5 | none in goldens; numeric options must be complete finite numbers (1e-3 and +0.1 now accepted); `--maf` in [0, 0.5] |
| 4955607 | U6 | none in goldens; one input format, matching reference, map required unless only nSL/XP-nSL; warnings for ignored flags |
| 903a394 | U12 | none in goldens; `--help` on stdout, `--version`, usage with no arguments, command line echoed once |
| 8d5d5c0 | notes | step 5 rows in this file |
| a1b775f | P4 | none; EHH kernels free their work arrays on every return path (p4 iHS leaked 37 MB; peak RSS 44.7 -> 9.1 MiB) |
| 067e9ab | P2 | none; phased data share one all-zero bitset (phased load RSS -31%) |
| 8a17e3a | P1 | none; VCF read through a 1 MiB gzread buffer (load about 2.7x faster) |
| cfa6f90 | B50 | VCF iHS/XP-EHH use `--map` also when nSL/XP-nSL is co-requested (was physical distance); goldens unchanged |
| b38368b | P1 | none; VCF positions taken from the haplotype pass (big50k load 1.09-1.29 -> 0.81-0.82 s) |
| 436e402 | P3 | none; gene annotation resumes its scan (1M SNPs, 2000 genes: 2.4 -> 2.1 s; 41.9 s on 3.0.1, mostly B48) |
| 3c79b85 | P6, B51 | none on well-formed input; norm read/normalise passes buffered (5M rows 10.8 -> 6.6 s); `--ihs-detail` column count from the first data line |
| 4c51bce | P5 | none; group splitting without per-step allocation (nSL 85.8 -> 14.9 s, single thread) |
| 47776f7 | B40 | none at the suite's sizes; exact pair counts (nCk(n,2) was off above 59,092 haplotypes; the factln cache race is gone) |

Each commit message gives the measured differences in full.

Correction to the 82d0042 message: its sentence "loaded 811 instead of 998
loci and reported 13 instead of 178 cores" pairs two different
comparisons. On the allele-swapped p1 at `--maf 0.1`, the pre-B32 build
loads 811 loci and reports 13 cores, and the fixed build loads 998 loci and
reports 198 cores. 178 is the pre-B32 build's count on the original,
unswapped p1, so the asymmetry was 178 cores unswapped against 13 swapped.

Corrections to later commit messages:
- f97de30 (B12) and 3558905 (B13) give their `make check` counts as from
  a "(clean build)"; no `make clean` preceded those runs, so they were
  incremental builds (header dependencies tracked since 341bdc0). At
  eea08a5, `make clean` followed by `make check` gives 155 passed, 3
  known failures, 0 failed.
- d2bd8ed (B24) calls the release-build out-of-bounds write hidden; it
  ran silently on a018ca0 but crashed (SIGTRAP) on a 22769f4 build.
- 47776f7 is titled "B52", but it fixes review finding B40 (the factln
  cache race, the only ThreadSanitizer report, and nCk's inexactness
  above 59,092); there is no B52. The pair counts no longer call nCk, so
  the race is gone by construction. TSan does not run on the macOS build
  machine, so no TSan rerun was made.

## Validation against an independent reference

The review wrote a Python reference implementation from the manual's
equations and the v2.0.3 loop semantics. Where a PI decision changed a
definition, a workspace copy of the reference was adjusted to that decision
(noted below). Printed precision is 6 significant digits. Runs on p1/p2
remove the 12 duplicated positions.

- Phased XP-EHH after B2: 2156/2156 sites, max abs difference 5.0e-7.
  3.0.1 misses the 14 restored sites and differs at 725. (Run after commit
  8e701b8, whose message says the check was still pending.)
- Unphased XP-EHH after B6 + B10: 2193/2193 sites, max abs 4.3e-6, with the
  reference's skip rule set to "one genotype class holds everyone".
- Unphased iHS after B11: 1419/1419 cores; ihh1, ihh0 and ihs within 5.0e-6.
- Phased nSL after B7 (no cutoff): example.hap with example.map 5374/5374,
  with the 300 kb gap map 5174/5174, within 5e-6 relative.
- iHS and nSL across the 100 kb gap: 4782/4782 and 5374/5374 (already
  correct before B8).
- iHH12 after B8: 5477/5477 within 4.6e-6 with the scale applied once. 3.0.1
  matched only when the reference squared the scale.
- iHH12 after Q5 (EHH12 truncation): example, 300 kb gap and 100 kb gap maps,
  each with and without `--trunc-ok`; all six site sets identical, values
  within 4.1e-6.
- iHS `--alt` after B9 (squared frequencies, stop when all haplotypes are
  unique): 998/998 cores within 4.8e-6.
- iHH12 `--alt`: 5205/5205 (5574 with `--trunc-ok`) within 4.6e-6, against
  the reference extended to squared-frequency EHH12.
- EHH12 output (`--ehh12`, binomial and `--alt`, derived allele major and
  minor at the core): every cell of 302, 275 and 934 rows matches an
  independent computation from haplotype counts.
- XP-nSL after B47, on p1/p2 restricted to selscan's retained set
  (duplicated positions and jointly monomorphic sites removed, 2229 sites):
  phased 2156/2156, unphased 2193/2193, same site sets; sL within 4.7e-5
  (rounding on values near 20), XP-nSL within 3.8e-6. The pre-B47 build
  differs by more than 1e-4 at 1427 phased and 2188 unphased sites.
- norm frequency bins after B34, against the normalisation track's
  refnorm.py (exact left-closed rule, 20 bins): p1 iHS 503/503, p1 nSL
  816/816, example.hap iHS 4782/4782, same locus sets, within 5e-6.
- `--crit-percent` after B19/U10, against a numpy computation: p4 iHS at
  100 bins 2.24509 / -2.32029 from 1419 loci (equal); p1 iHH12 one-tailed
  1.20591, 58 of 1145 loci flagged at 0.05.
- Window percentile labels after B22/U4, against refwin.py with the `>`
  rule: tie-heavy and random synthetic sets 0 mismatches of 1001 and
  3001; example.hap iHS and p1/p2 XP-EHH (both tails) 0 of 50.
- Gene table `*_adj` after B26: within 5.2e-6 of numpy least squares
  (30-gene synthetic XP table).
- Gene-set permutation test after B38: 0.306 (seed 7) against a
  200,000-permutation numpy estimate of 0.299. With the built-in
  generator (no GSL, step 7): 0.304 at seed 7, mean 0.2997 over seeds
  100-119 (sd 0.0043; the binomial sd at 10,000 permutations is 0.0046).
- UBSan with libc++ hardening (`_LIBCPP_HARDENING_MODE_DEBUG`) on
  8150dd1: the 38 norm checks pass with no reports.
- The same build on 3c79b85: the 39 norm checks pass. On the tree
  committed as 4c51bce (P5): all 183 checks pass.
- P5 and B40 (47776f7), outputs byte-identical to 436e402: phased iHS and nSL
  (1000 samples x 50k sites; p4.2000), unphased iHS, iHH12, XP-EHH and
  XP-nSL (phased and unphased, p1/p2), EHH and EHH12 at 5 loci each.
- P6: normalised output and log byte-identical to 436e402 on 5M-row iHS,
  1M-row XP-EHH and iHH12 tables, and on an `--ihs-detail` table with
  15-17-digit fields and nan/inf/-0/NaN/1.5x tokens.

v2.1.3 is not a baseline for the XP fixes. It drops the same 14
joint-frequency-0.5 sites as 3.0.1 and agrees with 3.0.1 exactly on every
site they share.

## Deliberate behaviour decisions (PI)

- B4: XP-EHH/XP-nSL/iHH12 keep p1 = 1.0 sites; they are completed sweeps in
  population 1 and are normalised in the single genome-wide bin. iHS/nSL
  exclude derived frequency 1.0, which has no frequency class, and norm
  reports the count.
- B10: the all-heterozygous case is also treated as "no variation", which
  goes beyond the review's suggestion and matches the unphased iHS kernel.
- U12: a missing `--ehh12` locus is an error, as for `--ehh`.
- B7: nSL uses no EHH cutoff; it is governed by `--max-extend-nsl`.
- Q5: iHH12 truncates when EHH12, not whole-sample EHH, falls below
  `--cutoff`. This lowers the number of reported cores (more reach a
  chromosome end; `--trunc-ok` keeps them) and raises every iHH12 value.
- B9: `--alt` is the sum of squared haplotype frequencies; integration
  stops when all haplotypes are unique.
- B32: a site exactly at `--maf` is kept, whichever allele is minor.
- B17/Q1: each statistic in a combined run equals its single run; XP is never
  MAF-filtered because iHS is co-requested.
- Duplicated positions: every reader drops all records that share a
  position and logs how many.
- Q7/B15: XP input is intersected by position, and REF/ALT must match;
  each kind of dropped record is counted in the log.
- B47 (found during B15, not in the review report): XP-nSL counts distance
  in retained sites, as nSL does. A definition change made for
  consistency; confirmed by the PI.
- Q9/B34: frequency bins are [b/B, (b+1)/B), computed as
  floor(f*B + 1e-6*B); the tolerance covers 6-digit printing.
- Q10/B19, U10, U4: `--crit-percent` over every locus with a finite
  normalised score (including bins below 20 loci, which are not
  written); iHH12 one-tailed; window labels use `>`; a warning when a
  SNP-count bin has fewer than 100 windows. Rank-based empirical
  percentiles remain a possible future change.
- Q11/B20: `--log-input` takes the bin count from the log; an explicit
  `--bins` must agree.
- Q13/U12: help on stdout with exit 0, `--version`, a missing `--ehh12`
  locus is an error.
- Step 5 choices made here, not by the PI: `--pmap` together with `--map`
  and `--wagh` without XP only warn; a string option followed by another
  option is an error; `--maf` outside [0, 0.5] is an error; `--multi-param`
  rejects unknown keys.
- Implementation choices made here, not by the PI: B25 writes NA in the
  gene table's `*_adj` columns when fewer than 3 genes pass `--min-snps`;
  B38 adds `--perm-seed` and reports (count+1)/(n+1); B18 treats any
  non-finite score or frequency as missing and logs the count; the norm
  log records bin means and variances at 17 significant digits.
- Step 6 findings and choices made here, not by the PI: B50 and B51
  were found during step 6 and are not in the review report. P5 visits
  groups in first-touched order; results do not depend on it (integer
  sums). P6 reports write errors on norm output (exit 1); they were
  ignored.

## Performance (step 6)

Wall time, release build, M-series Mac (8 cores), 7a93ff1 -> 47776f7.
Synthetic big50k: 1000 samples x 50k sites, 28 MiB VCF. p4.2000: 2000
sites x 3202 samples. Single runs, so +-5%. These compare what each
version computes, and the definitions differ: nSL now integrates without
an EHH cutoff (B7) and iHH12 truncates on EHH12 (Q5), so both do more
work than in 3.0.1.

| run | 7a93ff1 | 47776f7 | ratio |
|---|---|---|---|
| iHS big50k t1 | 9.35 s | 2.77 s | 3.4x |
| iHS big50k t8 | 4.43 s | 1.13 s | 3.9x |
| nSL big50k t8 | 21.15 s | 2.92 s | 7.2x |
| iHS p4 t1 | 21.32 s | 4.73 s | 4.5x |
| iHH12 p4 t1 | 5.82 s | 4.21 s | 1.4x |
| iHS p4 unphased t1 | 14.7 s | 3.69 s | 4.0x |
| norm iHS 5M rows | 16.2 s | 6.66 s | 2.4x |

Per-commit measurements are in the commit messages. Main contributions:
VCF loading (P1, 3.1-3.4 -> 0.81-0.82 s on big50k), group splitting
(P5), exact pair counts (B40), norm I/O (P6). Peak memory: P4 and P2 (p4
iHS, 4 threads, 44.7 -> 7.0 MiB).

Tried and not committed: writing rows with "\n" instead of `endl` (one
flush per row). No measurable change on local disk; the writes overlap
the computation. It may matter on network filesystems.

## Known open issues

- KNOWN_FAIL in tests/run_tests.sh is empty since 02aaea8 (B39).
- B43 (phased XP-EHH `--alt`) differs from the reference at 862 of 2083
  sites. Part of it is explained: XP integrates each population only while
  that population is not all-unique, and under `--alt` the all-unique EHH
  floor is 1/n, not 0. That accounts for 53 sites; 809 are unexplained.
  Binomial XP-EHH matches exactly.
- iHH12 `--alt` code mixes squared and binomial terms for the top two
  haplotypes. Left alone deliberately: the EHH12 adjustment is 2*t1*t2 under
  either measure, and iHH12 `--alt` was identical before and after the
  equivalent change.
- The loader differences found while building the harness (duplicated
  positions, jointly monomorphic XP sites) are resolved by f985491 and
  28f0f22; tped-equiv and thap-equiv now run on the unfiltered p1/p2.
- `--filter-vcf` (the multi-allelic filter) tests REF/ALT for single
  A/C/G/T bases, so on 0/1-coded VCFs such as example/ it removes every
  record (B39, filter-vcf-01).
- No golden case runs unphased `--wagh`, so that part of feec363 is
  untested.
- B24 is a regression check only on release builds: the out-of-bounds
  write ran silently on a018ca0 but crashed (SIGTRAP) on a 22769f4
  build, and a libc++-hardened build aborts on it.
- norm items from the review not yet addressed: B35 (`--min-snps 0`), B36
  (window analysis assumes one sorted chromosome per file), B46
  (`--win-offset` counts leading SNPs), U3 (`--ihs-detail` input), U11
  (empty inputs), M11 (help and manual).
- norm speed after P6: about 2.5 of 6.6 s on 5M rows is `%.6g`
  formatting. `std::to_chars` would remove most of it but needs macOS
  13.3 or later in Apple's libc++; this waits for the release deployment
  target (step 7). The separate counting pass (checkIHSfile and its XP
  and iHH12 versions, about 0.9 s) remains.
- binom.cpp (nCk, factln) is no longer used by the kernels since 47776f7 (B40).
- iHS/nSL gene-table rows end with an extra tab (one more field than the
  header); left as is.
- Makefile_win and the unused Makefile_macos-arm64/-x86-64 still lack
  header dependencies (M4).
- `--multi-chr` is deferred until after the bug fixes; it is a new feature
  to add.
