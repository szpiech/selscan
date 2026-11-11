# Branch `code-review`: review fixes on top of 3.0.1

Base: `main` at 7a93ff1 ("Patch 3_0_1"). Finding ids (B1, U12, ...) refer to
the code review report (selscan_code_review.md). Every commit builds from a
clean tree and passes `make check`. Run `make clean` before testing a change
to a header, because the Makefiles have no header dependencies (M1, not yet
fixed).

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

Each commit message gives the measured differences in full.

## Validation against an independent reference

The review wrote a Python reference implementation from the manual's
equations and the v2.0.3 loop semantics. The runs below use example p1/p2
with the 12 duplicated positions removed; printed precision is 6 significant
digits.

- Phased XP-EHH after B2: 2156/2156 sites, max abs difference 5.0e-7.
  3.0.1 misses the 14 restored sites and differs at 725. (Run after commit
  8e701b8, whose message says the check was still pending.)
- Unphased XP-EHH after B6 + B10: 2193/2193 sites, max abs 4.3e-6, with the
  reference's skip rule set to "one genotype class holds everyone".
- Unphased iHS after B11: 1419/1419 cores; ihh1, ihh0 and ihs within 5.0e-6.

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

## Known open issues

- Still failing on purpose, listed in tests/run_tests.sh KNOWN_FAIL: B27
  (`--keep-low-freq` ignored in combined runs), B32 (MAF boundary asymmetry,
  `--multi-param` differs from single runs) and B39 (uncaught exceptions
  abort with 134).
- Found while building the harness, not in the review report: the VCF loader
  drops sites at duplicated positions and the TPED loader keeps them; in XP
  mode only the VCF loader drops sites monomorphic in the two panels
  combined. The extra sites change iHH through the integration steps and gap
  scaling, so the same data gives different XP output from TPED and VCF.
  The tped-equiv fixture removes both kinds of site so that it tests file
  reading alone.
- No golden case runs unphased `--wagh`, so that part of feec363 is
  untested.
