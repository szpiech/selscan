#!/bin/bash
# selscan regression suite.
#
#   ./tests/run_tests.sh                 # check src/selscan against tests/expected/
#   ./tests/run_tests.sh --regen         # re-record tests/expected/ from the current build
#   ./tests/run_tests.sh --bin PATH      # test another binary
#   ./tests/run_tests.sh --rtol 1e-5     # loosen the float tolerance (default 1e-6)
#   ./tests/run_tests.sh --keep          # keep the scratch directory
#   ./tests/run_tests.sh --list          # list case names
#   ./tests/run_tests.sh <case> [...]    # run only the named cases
#   (or `make check` from src/)
#
# Three kinds of case:
#   golden     a statistic's output table compared column by column with the
#              recorded table in tests/expected/ (header and integer columns
#              exactly, floats to --rtol). A golden records what the code DOES;
#              a commit that changes one must say so and show the difference.
#   invariant  a property that needs no golden: thread count does not change
#              output, a combined run equals separate runs, XP statistics are
#              antisymmetric under swapping the populations, TPED input equals
#              the same data as VCF, and so on. These encode what the code
#              SHOULD do, and are where the review findings show up.
#   status     bad input must exit 1 with a message, never crash.
#
# Known failures: KNOWN_FAIL lists checks that fail on purpose, each tagged
# with the finding in selscan_code_review.md that explains it. They print as
# "xfail" and do not fail the run. A listed check that PASSES is an error
# ("XPASS"), so a fix must remove its own entries in the same commit.
#
# All fixtures are derived at run time from files already in example/; no
# test data is committed beyond the expected tables.

set -u

TESTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$TESTS_DIR/.." && pwd)"
EXAMPLE="$ROOT_DIR/example"
BIN="$ROOT_DIR/src/selscan"
EXPECTED="$TESTS_DIR/expected"
REGEN=0; KEEP=0; LIST=0; RTOL=1e-6
SELECTED=()

while [ $# -gt 0 ]; do
    case "$1" in
        --regen) REGEN=1; shift ;;
        --keep)  KEEP=1; shift ;;
        --list)  LIST=1; shift ;;
        --bin)   BIN="$2"; shift 2 ;;
        --rtol)  RTOL="$2"; shift 2 ;;
        --help|-h) sed -n '2,32p' "$BASH_SOURCE"; exit 0 ;;
        -*) echo "unknown option: $1" >&2; exit 64 ;;
        *)  SELECTED+=("$1"); shift ;;
    esac
done
case "$BIN" in /*) ;; *) BIN="$PWD/$BIN" ;; esac

# ------------------------------------------------------------ known failures
# <check name> <finding id>. Remove an entry in the commit that fixes it.
KNOWN_FAIL="
maxgap-effective:ihs  B1
maxgap-effective:nsl  B1
xp-fixed-diff:xpehh  B2
xp-fixed-diff:xpehh-unphased  B2
norm-freq1-excluded:ihs  B4
norm-single-bin-linear:xpehh  B4
norm-single-bin-linear:xpnsl  B4
xp-swap-antisymmetry:xpehh-unphased  B6
xp-swap-antisymmetry:xpnsl-unphased  B6
multistat==single:ihs--keep-low-freq  B27
multistat==single:nsl--keep-low-freq  B27
multistat==single:ihh12--keep-low-freq  B27
multiparam==single:pconfig2  B32
multiparam==single:pconfig7  B32
multiparam==single:pconfig9  B32
clean-error:all-filtered  B39
clean-error:filter-vcf-01  B39
clean-error:bad-json  B39
"
known_fail() { printf '%s\n' "$KNOWN_FAIL" | awk -v k="$1" '$1==k {print $2; exit}'; }

# ---------------------------------------------------------------- reporting
NPASS=0; NFAIL=0; NXFAIL=0
FAILED=()
if [ -t 1 ]; then G=$'\033[32m'; R=$'\033[31m'; Y=$'\033[33m'; N=$'\033[0m'; else G=; R=; Y=; N=; fi
pass() {
    local id; id="$(known_fail "$1")"
    if [ -n "$id" ] && [ "$REGEN" = 0 ]; then
        printf '  %sXPASS%s %s: listed as known failure %s but passes; remove it from KNOWN_FAIL\n' "$R" "$N" "$1" "$id"
        NFAIL=$((NFAIL+1)); FAILED+=("$1")
    else
        printf '  %sok%s    %s\n' "$G" "$N" "$1"; NPASS=$((NPASS+1))
    fi
}
fail() {
    local id; id="$(known_fail "$1")"
    if [ -n "$id" ]; then
        printf '  %sxfail%s %s [%s]: %s\n' "$Y" "$N" "$1" "$id" "$2"; NXFAIL=$((NXFAIL+1))
    else
        printf '  %sFAIL%s  %s: %s\n' "$R" "$N" "$1" "$2"; NFAIL=$((NFAIL+1)); FAILED+=("$1")
    fi
}

selected() {
    [ ${#SELECTED[@]} -eq 0 ] && return 0
    local c; for c in "${SELECTED[@]}"; do [ "$c" = "$1" ] && return 0; done
    return 1
}

# run <log> args... : run the binary, record its exit status in RC
RC=0
# (in a subshell, so a crash report from the shell lands in the log too)
run() { local log="$1"; shift; ( "$BIN" "$@"; exit $? ) > "$log" 2>&1; RC=$?; }

# compare_table <check> <actual> <golden-name> [noheader]
# Header must match exactly, integer-looking fields exactly, floats to RTOL.
# The failure message names the column, the number of rows and the worst
# relative deviation.
compare_table() {
    local name="$1" actual="$2" gname="$3" nohdr="${4:-}"
    local golden="$EXPECTED/$gname.tsv.gz"
    if [ ! -s "$actual" ]; then fail "$name" "no output at $(basename "$actual")"; return; fi
    if [ "$REGEN" = 1 ]; then gzip -n -c "$actual" > "$golden"; pass "$name (recorded $gname.tsv.gz)"; return; fi
    if [ ! -f "$golden" ]; then fail "$name" "no golden $gname.tsv.gz (run --regen)"; return; fi
    local report
    report=$(awk -v rtol="$RTOL" -v nohdr="$nohdr" '
        function isnum(x) { return (x ~ /^[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?$/) }
        function isint(x) { return (x ~ /^[+-]?[0-9]+$/) }
        FNR==NR {
            if (FNR==1 && nohdr=="") { ehdr=$0; en=NF; for(i=1;i<=NF;i++) ename[i]=$i; next }
            if (FNR==1) { en=NF; for(i=1;i<=NF;i++) ename[i]="col" i }
            for(i=1;i<=NF;i++) E[FNR,i]=$i; erows=FNR; enf[FNR]=NF
            next
        }
        {
            if (FNR==1 && nohdr=="") { if ($0 != ehdr) { print "header differs"; bad=1; exit } next }
            if (FNR > erows) { print "actual has more rows than golden"; bad=1; exit }
            if (NF != enf[FNR]) { print "row " FNR ": " NF " fields, golden has " enf[FNR]; bad=1; exit }
            for (i=1;i<=NF;i++) {
                e=E[FNR,i]; a=$i
                if (e==a) continue
                if (!isnum(e) || !isnum(a) || (isint(e) && isint(a))) { diffcol[i]++; if (!(i in first)) first[i]=FNR; continue }
                d = (e==0) ? (a==0?0:1) : (e-a)/e; if (d<0) d=-d
                if (d > maxd[i]) maxd[i]=d
                if (d > rtol) { diffcol[i]++; if (!(i in first)) first[i]=FNR }
            }
            arows=FNR
        }
        END {
            if (bad) exit
            if (arows != erows) { print "row count differs: golden " erows " actual " arows; exit }
            msg=""
            for (i=1;i<=en;i++) if (diffcol[i])
                msg = msg sprintf("%s: %d rows differ, first at line %d (max rel %.3g); ", ename[i], diffcol[i], first[i], maxd[i])
            if (msg != "") print msg
        }
    ' <(gzip -dc "$golden") "$actual")
    if [ -z "$report" ]; then pass "$name"; else fail "$name" "$report"; fi
}

# identical <check> <a> <b> : byte-identical outputs
identical() {
    if [ ! -s "$2" ] || [ ! -s "$3" ]; then fail "$1" "missing output"; return; fi
    if cmp -s "$2" "$3"; then pass "$1"; else fail "$1" "outputs differ ($(diff "$2" "$3" | grep -c '^[<>]') lines)"; fi
}

# ---------------------------------------------------------------- case list
CASES="ihs ihs-detail ihs-alt ihs-keep nsl xpehh xpehh-wagh xpnsl ihh12 pi ehh ehh12
ihs-unphased nsl-unphased xpehh-unphased xpnsl-unphased ehh-unphased
ihs-hap nsl-thap xpehh-hap ihs-tped ihs-gap nsl-gap ihh12-gap
norm-ihs norm-nsl norm-xpehh norm-ihh12 norm-ihs-bpwin
threads multistat multiparam xp-swap maxgap xp-fixed-diff tped-equiv norm-invariants errors help"

if [ "$LIST" = 1 ]; then for c in $CASES; do echo "$c"; done; exit 0; fi
if [ ! -x "$BIN" ]; then echo "ERROR: no binary at $BIN (run make in src/ first)." >&2; exit 70; fi
for f in p1.vcf.gz p2.vcf.gz example.hap example.thap example.map example2.p1.hap example2.p2.hap example2.map p1.tped param2.json; do
    [ -f "$EXAMPLE/$f" ] || { echo "ERROR: missing example/$f" >&2; exit 70; }
done

WORK="$(mktemp -d "${TMPDIR:-/tmp}/selscan-tests.XXXXXX")"
cleanup() { if [ "$KEEP" = 1 ]; then echo "scratch kept at $WORK"; else rm -rf "$WORK"; fi; }
trap cleanup EXIT
mkdir -p "$EXPECTED"
cd "$WORK" || exit 70
for f in "$EXAMPLE"/*; do ln -s "$f" .; done

# ---------------------------------------------------------------- fixtures
# gap.map: example.map with a 300 kb gap (and 0.3 cM) inserted before the
# middle locus, larger than the default --max-gap of 200 kb.
awk '{ if (NR > 4162) { $4 += 300000; $3 = sprintf("%.15g", $3 + 0.3) } print }' example.map > gap.map

# fd1/fd2.vcf: p1/p2 with every 400th site (from the 200th) turned into a fixed
# difference, all 0|0 in p1 and all 1|1 in p2. With equal panel sizes these
# sites sit at joint frequency exactly 0.5 and carry no heterozygotes.
for k in 1 2; do
    gt="0|0"; [ $k = 2 ] && gt="1|1"
    gzip -dc p$k.vcf.gz | awk -v gt="$gt" 'BEGIN{OFS="\t"} /^#/ {print; next}
        { r++; if (r % 400 == 200) { for (i=10;i<=NF;i++) $i=gt; print $3 > "fd_ids.txt" } print }' > fd$k.vcf
done

# dd1/dd2.vcf: p1/p2 restricted to records that both loaders keep, and
# c1/c2.tped: the same records converted to TPED, to check that TPED input is
# read exactly like the same data in VCF. Two filters differ between the
# loaders and are removed here so this case tests reading alone: the VCF
# loader drops sites at duplicated positions (12 in p1) and, in XP mode,
# sites monomorphic in the two panels combined; the TPED loader does neither.
# Sites at joint frequency 0.5 are also removed (dropped by VCF XP, B2).
gzip -dc p1.vcf.gz | awk '!/^#/ {print $2}' | uniq -d > dup_pos.txt
for k in 1 2; do
    gzip -dc p$k.vcf.gz | awk '!/^#/ { c=0; for (i=10;i<=NF;i++) c += (substr($i,1,1)=="1") + (substr($i,3,1)=="1"); print c, 2*(NF-9) }' > cnt$k.txt
done
paste -d' ' cnt1.txt cnt2.txt | awk '{ t=$1+$3; n=$2+$4; print (t==0 || t==n || 2*t==n) ? 1 : 0 }' > drop_row.txt
for k in 1 2; do
    gzip -dc p$k.vcf.gz | awk 'FILENAME==ARGV[1] {d[$1]=1; next} FILENAME==ARGV[2] {x[FNR]=$1; next}
        /^#/ {print; next} { r++; if (!($2 in d) && !x[r]) print }' dup_pos.txt drop_row.txt - > dd$k.vcf
    awk '/^#/ {next} { printf "%s %s %.10g %s", $1, $3, $2/1e6, $2;
        for (i=10;i<=NF;i++) printf " %s %s", substr($i,1,1), substr($i,3,1); printf "\n" }' dd$k.vcf > c$k.tped
done

echo "selscan regression suite"
echo "  binary : $BIN"
echo "  rtol   : $RTOL"
[ "$REGEN" = 1 ] && echo "  mode   : RECORDING goldens into $EXPECTED"
echo

# ---------------------------------------------------------------- golden cases
# golden <case> <output-suffix> [noheader] -- args...
golden() {
    local name="$1" out="$2" nohdr=""; shift 2
    if [ "$1" = noheader ]; then nohdr=1; shift; fi
    [ "$1" = -- ] && shift
    selected "$name" || return 0
    run "$name.log" "$@" --out "g_$name"
    # A crash after the table is written still leaves a comparable table;
    # the exit status is checked separately so the two failures are distinct.
    if [ "$RC" -ne 0 ]; then fail "status:$name" "exit $RC"; else pass "status:$name"; fi
    compare_table "golden:$name" "g_$name.$out" "$name" "$nohdr"
}

golden ihs            ihs.out          -- --ihs --vcf p1.vcf.gz --pmap
golden ihs-detail     ihs.out          -- --ihs --ihs-detail --vcf p1.vcf.gz --pmap
golden ihs-alt        ihs.alt.out      -- --ihs --alt --trunc-ok --vcf p1.vcf.gz --pmap
golden ihs-keep       ihs.out          -- --ihs --keep-low-freq --vcf p1.vcf.gz --pmap
golden nsl            nsl.out          -- --nsl --vcf p1.vcf.gz
golden xpehh          xpehh.out        -- --xpehh --vcf p1.vcf.gz --vcf-ref p2.vcf.gz --pmap
golden xpehh-wagh     xpehh.out        -- --xpehh --wagh --vcf p1.vcf.gz --vcf-ref p2.vcf.gz --pmap
golden xpnsl          xpnsl.out        -- --xpnsl --vcf p1.vcf.gz --vcf-ref p2.vcf.gz
golden ihh12          ihh12.out        -- --ihh12 --vcf p1.vcf.gz --pmap
golden pi             pi.100bp.out noheader -- --pi --vcf p1.vcf.gz --pmap
golden ehh            ehh.chr1:29243.out   -- --ehh chr1:29243 --vcf p1.vcf.gz --pmap
golden ehh12          ehh12.chr1:29243.out -- --ehh12 chr1:29243 --vcf p1.vcf.gz --pmap
golden ihs-unphased   ihs.out          -- --ihs --unphased --vcf p1.vcf.gz --pmap
golden nsl-unphased   nsl.out          -- --nsl --unphased --vcf p1.vcf.gz
golden xpehh-unphased xpehh.out        -- --xpehh --unphased --vcf p1.vcf.gz --vcf-ref p2.vcf.gz --pmap
golden xpnsl-unphased xpnsl.out        -- --xpnsl --unphased --vcf p1.vcf.gz --vcf-ref p2.vcf.gz
golden ehh-unphased   ehh.chr1:29243.out   -- --ehh chr1:29243 --unphased --vcf p1.vcf.gz --pmap
golden ihs-hap        ihs.out          -- --ihs --hap example.hap --map example.map
golden nsl-thap       nsl.out          -- --nsl --thap example.thap --map example.map
golden xpehh-hap      xpehh.out        -- --xpehh --hap example2.p1.hap --ref example2.p2.hap --map example2.map
golden ihs-tped       ihs.out          -- --ihs --tped p1.tped
golden ihs-gap        ihs.out          -- --ihs --hap example.hap --map gap.map
golden nsl-gap        nsl.out          -- --nsl --hap example.hap --map gap.map
golden ihh12-gap      ihh12.out        -- --ihh12 --hap example.hap --map gap.map

# norm cases read the RECORDED golden statistic tables, so they test norm alone.
norm_golden() {
    local name="$1" src="$2" out="$3"; shift 3
    selected "$name" || return 0
    if [ ! -f "$EXPECTED/$src.tsv.gz" ]; then fail "golden:$name" "needs golden $src"; return; fi
    gzip -dc "$EXPECTED/$src.tsv.gz" > "n_$name.in"
    run "$name.log" norm "$@" --files "n_$name.in"
    if [ "$RC" -ne 0 ]; then fail "status:$name" "exit $RC"; else pass "status:$name"; fi
    compare_table "golden:$name" "n_$name.in.$out" "$name"
}
norm_golden norm-ihs       ihs      20bins.norm --ihs --bins 20
norm_golden norm-nsl       nsl      20bins.norm --nsl --bins 20
norm_golden norm-xpehh     xpehh    norm        --xpehh
norm_golden norm-ihh12     ihh12    norm        --ihh12
norm_golden norm-ihs-bpwin ihs-hap  20bins.norm --ihs --bins 20 --bp-win --winsize 100000
if selected norm-ihs-bpwin; then
    compare_table "golden:norm-ihs-bpwin-windows" n_norm-ihs-bpwin.in.20bins.norm.100kb.windows norm-ihs-bpwin-windows noheader
fi

[ "$REGEN" = 1 ] && { echo; echo "goldens recorded in $EXPECTED"; exit 0; }

# ---------------------------------------------------------------- invariants
if selected threads; then
    for st in "ihs|--ihs --vcf p1.vcf.gz --pmap" "nsl|--nsl --vcf p1.vcf.gz" \
              "xpehh|--xpehh --vcf p1.vcf.gz --vcf-ref p2.vcf.gz --pmap" \
              "xpnsl|--xpnsl --unphased --vcf p1.vcf.gz --vcf-ref p2.vcf.gz" \
              "ihh12|--ihh12 --vcf p1.vcf.gz --pmap" "ihs-unphased|--ihs --unphased --vcf p1.vcf.gz --pmap"; do
        n=${st%%|*}; a=${st#*|}; e=${n%%-*}
        run t1.log $a --threads 1 --out t1_$n; run t4.log $a --threads 4 --out t4_$n
        identical "threads1==threads4:$n" t1_$n.$e.out t4_$n.$e.out
    done
fi

if selected multistat; then
    # Every statistic in a combined run must equal the same statistic run alone.
    for keep in "" "--keep-low-freq"; do
        run m.log --ihs --nsl --ihh12 --vcf p1.vcf.gz --pmap $keep --out multi
        for e in ihs nsl ihh12; do
            run s.log --$e --vcf p1.vcf.gz --pmap $keep --out single
            identical "multistat==single:$e$keep" multi.$e.out single.$e.out
        done
    done
fi

if selected multiparam; then
    # Each parameter set of a --multi-param run must equal the equivalent single run.
    run mp.log --ihs --vcf p1.vcf.gz --pmap --multi-param param2.json --out mp
    awk '
        /"params"/ { inp=1 }
        inp && /\{/ && !/"params"/ { i++; maf[i]=0.05; cut[i]=0.05; me[i]=1000000 }
        inp && /"maf"/        { v=$0; gsub(/.*: *|[ ,}]/,"",v); maf[i]=v }
        inp && /"cutoff"/     { v=$0; gsub(/.*: *|[ ,}]/,"",v); cut[i]=v }
        inp && /"max-extend"/ { v=$0; gsub(/.*: *|[ ,}]/,"",v); me[i]=v }
        END { for (j=1;j<=i;j++) print j-1, maf[j], cut[j], me[j] }' param2.json > mp_cases.txt
    while read -r i maf cut me; do
        run sp.log --ihs --vcf p1.vcf.gz --pmap --maf "$maf" --cutoff "$cut" --max-extend "$me" --out sp
        identical "multiparam==single:pconfig$i" mp.pconfig$i.ihs.out sp.ihs.out
    done < mp_cases.txt
fi

# xp_antisym <check> <ab> <ba> : swapping query and reference negates the
# statistic and covers the same sites.
xp_antisym() {
    local r
    r=$(awk -F'\t' 'FNR==1 {next} FNR==NR {a[$2]=$NF; na++; next}
        { nb++; if (!($2 in a)) only++; else { d=a[$2]+$NF; if (d<0) d=-d; if (d>w) w=d } }
        END { if (na != nb || only) printf "site sets differ (%d vs %d); ", na, nb; if (w > 1e-5) printf "max |AB+BA| = %.4g", w }' "$2" "$3")
    if [ -z "$r" ] && [ -s "$2" ]; then pass "$1"; else fail "$1" "${r:-no output}"; fi
}
if selected xp-swap; then
    for ph in "" "--unphased"; do for st in xpehh xpnsl; do
        pm=; [ $st = xpehh ] && pm=--pmap
        run ab.log --$st $ph --vcf p1.vcf.gz --vcf-ref p2.vcf.gz $pm --out ab
        run ba.log --$st $ph --vcf p2.vcf.gz --vcf-ref p1.vcf.gz $pm --out ba
        xp_antisym "xp-swap-antisymmetry:$st${ph:+-unphased}" ab.$st.out ba.$st.out
    done; done
fi

if selected maxgap; then
    # gap.map has a 300 kb gap: with the default --max-gap 200000 every core
    # whose haplotypes reach the gap must be dropped, so the output must have
    # fewer rows than a run with the gap limit effectively disabled.
    for st in ihs nsl; do
        run g1.log --$st --hap example.hap --map gap.map --out mg_def
        run g2.log --$st --hap example.hap --map gap.map --max-gap 1000000000 --out mg_off
        a=$(($(wc -l < mg_def.$st.out))); b=$(($(wc -l < mg_off.$st.out)))
        if [ "$a" -lt "$b" ]; then pass "maxgap-effective:$st"
        else fail "maxgap-effective:$st" "$a rows with --max-gap 200000, $b with the gap limit disabled"; fi
    done
fi

if selected xp-fixed-diff; then
    # Fixed differences between the populations are the strongest XP signal and
    # must never be filtered as monomorphic.
    nfd=$(wc -l < fd_ids.txt | tr -d ' ')
    for ph in "" "--unphased"; do
        run fd.log --xpehh $ph --vcf fd1.vcf --vcf-ref fd2.vcf --pmap --out fd
        got=$(awk -F'\t' 'FNR==NR {w[$1]=1; next} FNR>1 && ($2 in w) {n++} END {print n+0}' fd_ids.txt fd.xpehh.out)
        if [ "$got" = "$nfd" ]; then pass "xp-fixed-diff:xpehh${ph:+-unphased}"
        else fail "xp-fixed-diff:xpehh${ph:+-unphased}" "$got of $nfd fixed-difference sites reported"; fi
    done
fi

if selected tped-equiv; then
    # The same data as TPED and as VCF must give identical results.
    run te1.log --ihs --tped c1.tped --pmap --out te_t; run te2.log --ihs --vcf dd1.vcf --pmap --out te_v
    identical "tped-equiv:ihs" te_t.ihs.out te_v.ihs.out
    run te3.log --xpehh --tped c1.tped --tped-ref c2.tped --pmap --out te_t
    run te4.log --xpehh --vcf dd1.vcf --vcf-ref dd2.vcf --pmap --out te_v
    identical "tped-equiv:xpehh" te_t.xpehh.out te_v.xpehh.out
fi

if selected norm-invariants; then
    # XP-EHH/XP-nSL are normalised in one genome-wide bin: every output row must
    # satisfy norm = (x - mean)/sd for the same mean and sd, and every input row
    # must be written. Checked by fitting the line through the smallest and
    # largest raw values and testing every row against it.
    for st in xpehh xpnsl; do
        gzip -dc "$EXPECTED/$st.tsv.gz" > nl_$st.out
        run nl.log norm --$st --files nl_$st.out
        r=$(awk -F'\t' -v st="$st" -v nin="$(($(wc -l < nl_$st.out)))" '
            NR==1 { for (i=1;i<=NF;i++) { if ($i==st) xc=i; if ($i=="norm_" st) nc=i }; next }
            { x[NR]=$xc+0; y[NR]=$nc+0; n=NR; if (lo=="" || x[NR]<x[lo]) lo=NR; if (hi=="" || x[NR]>x[hi]) hi=NR }
            END {
                if (n != nin) printf "%d of %d rows written; ", n-1, nin-1
                s=(y[hi]-y[lo])/(x[hi]-x[lo]); c=y[lo]-s*x[lo]
                for (i=2;i<=n;i++) { d=y[i]-(s*x[i]+c); if (d<0) d=-d; if (d > 1e-4*(1+(y[i]<0?-y[i]:y[i]))) bad++ }
                if (bad) printf "%d rows off the normalisation line", bad }' nl_$st.out.norm)
        if [ -z "$r" ]; then pass "norm-single-bin-linear:$st"; else fail "norm-single-bin-linear:$st" "$r"; fi
    done
    # iHS: loci at derived frequency 1.0 cannot be normalised in a frequency bin;
    # they must be left out of the output and reported, not written with another
    # locus's score.
    gzip -dc "$EXPECTED/ihs.tsv.gz" | awk 'BEGIN{OFS="\t"} NR>1 && NR%50==0 {$4="1"} {print}' > f1.ihs.out
    nf1=$(awk 'NR>1 && $4==1' f1.ihs.out | wc -l | tr -d ' ')
    run f1.log norm --ihs --bins 20 --files f1.ihs.out
    got=$(awk 'NR>1 && $4==1' f1.ihs.out.20bins.norm | wc -l | tr -d ' ')
    if [ "$got" = 0 ] && grep -qi "excluded $nf1 loci" f1.log; then pass "norm-freq1-excluded:ihs"
    else fail "norm-freq1-excluded:ihs" "$got of $nf1 freq=1 rows written; exclusion message $(grep -qi excluded f1.log && echo present || echo absent)"; fi
fi

if selected errors; then
    printf '{"params":[{"maf":0.1,}]}' > bad.json
    while IFS='|' read -r name args; do
        [ -z "$name" ] && continue
        run e.log $args --out err_$name
        if [ "$RC" -eq 1 ]; then pass "clean-error:$name"; else fail "clean-error:$name" "exit $RC (expected 1)"; fi
    done <<'ERRS'
no-stat|--vcf p1.vcf.gz --pmap
xp-no-ref|--xpehh --vcf p1.vcf.gz --pmap
no-map|--ihs --vcf p1.vcf.gz
all-filtered|--ihs --vcf p1.vcf.gz --pmap --maf 0.7
filter-vcf-01|--ihs --filter-vcf --vcf p1.vcf.gz --pmap
bad-locus|--ehh no_such_locus --vcf p1.vcf.gz --pmap
bad-locus-ehh12|--ehh12 no_such_locus --vcf p1.vcf.gz --pmap
bad-json|--ihs --vcf p1.vcf.gz --pmap --multi-param bad.json
ERRS
fi

if selected help; then
    run h.log --help
    if [ "$RC" -eq 0 ]; then pass "help-exit-0"; else fail "help-exit-0" "exit $RC"; fi
fi

# ---------------------------------------------------------------- summary
echo
printf 'passed %d, known failures %d, failed %d\n' "$NPASS" "$NXFAIL" "$NFAIL"
if [ "$NFAIL" -gt 0 ]; then printf 'failed: %s\n' "${FAILED[*]}"; exit 1; fi
exit 0
