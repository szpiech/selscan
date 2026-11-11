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

Each commit message gives the measured differences in full.

Correction to the 82d0042 message: its sentence "loaded 811 instead of 998
loci and reported 13 instead of 178 cores" pairs two different
comparisons. On the allele-swapped p1 at `--maf 0.1`, the pre-B32 build
loads 811 loci and reports 13 cores, and the fixed build loads 998 loci and
reports 198 cores. 178 is the pre-B32 build's count on the original,
unswapped p1, so the asymmetry was 178 cores unswapped against 13 swapped.

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
  in retained sites, as nSL does. This is a definition change made for
  consistency and awaits PI confirmation; it is a separate commit so it can
  be dropped.

## Known open issues

- Still failing on purpose, listed in tests/run_tests.sh KNOWN_FAIL: B39
  (uncaught exceptions abort with exit 134: all loci filtered, `--filter-vcf`
  on 0/1-coded VCF, malformed `--multi-param` JSON).
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
- Makefile_win and the unused Makefile_macos-arm64/-x86-64 still lack
  header dependencies (M4).
- `--multi-chr` is deferred until after the bug fixes; it is a new feature
  to add.
