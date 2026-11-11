# selscan regression suite

```
cd src && make check                         # build, then run the suite on ./selscan
./tests/run_tests.sh --bin path/to/selscan   # any binary
./tests/run_tests.sh ihs xp-swap             # selected cases (see --list)
./tests/run_tests.sh --regen                 # re-record tests/expected/ (see below)
```

Runtime is about 20 s on a laptop. Needs only bash, awk, gzip, cmp and diff.

## What is checked

- **golden**: 30 output tables (every statistic, phased and unphased; VCF,
  hap, thap and TPED input; `selscan norm` for each mode and for `--bp-win`)
  compared column by column with `tests/expected/*.tsv.gz`. Headers and
  integer fields must match exactly, floats to `--rtol` (default 1e-6). A
  failure names the column, the number of rows and the worst relative
  deviation. Each golden case also checks the exit status separately
  (`status:<case>`), so "crashed after writing correct output" and "wrong
  output" are different failures.
- **invariant**: properties that need no golden. Identical output at 1 and 4
  threads; every statistic in a combined run equal to the same statistic run
  alone; `--multi-param` sets equal to the equivalent single runs; XP-EHH/XP-nSL
  antisymmetric when the populations are swapped; `--max-gap` actually removes
  cores across a 300 kb gap; fixed differences between populations reach the XP
  output; the same data as TPED and as VCF give identical results; single-bin
  normalisation is linear and writes every row; frequency-1.0 loci are excluded
  from iHS normalisation and reported.
- **status**: bad input exits 1 with a message rather than crashing.

Fixtures (a map with an inserted gap, VCFs with planted fixed differences, TPED
converted from VCF) are derived at run time from `example/`. No test data is
committed apart from the expected tables.

## Known failures

`KNOWN_FAIL` at the top of `run_tests.sh` lists checks that fail on purpose,
each tagged with the finding in the code review that explains it. They print as
`xfail` and do not fail the run. A listed check that passes prints `XPASS` and
DOES fail the run, so the commit that fixes a finding must delete its entries.

## Changing a golden

The goldens were recorded from the unmodified v3.0.1 source (7a93ff1). A
commit that intentionally changes output re-records the affected goldens with
`--regen` and states in its message which tables changed and by how much. The
`norm-*` goldens are computed from the recorded statistic tables (not from a
fresh run), so a change to a statistic only moves its own golden, unless it
changes a table that norm reads; regenerate both in that case.
