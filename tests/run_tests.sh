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
CASES="ihs ihs-detail ihs-alt ihs-keep nsl xpehh xpehh-wagh xpnsl ihh12 pi ehh ehh12 ehh12-alt ehh12-minor
ihs-unphased nsl-unphased xpehh-unphased xpnsl-unphased ehh-unphased
ihs-hap nsl-thap xpehh-hap ihs-tped ihs-gap nsl-gap ihh12-gap ihs-gap100 nsl-gap100 ihh12-gap100
norm-ihs norm-nsl norm-xpehh norm-ihh12 norm-ihs-bpwin
threads multistat multiparam xp-swap maxgap xp-fixed-diff tped-equiv norm-invariants maf-symmetry pi-edges errors help"

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
# gap100.map: the same with a 100 kb gap (and 0.1 cM), between --gap-scale
# (20 kb) and --max-gap (200 kb), so cores cross it with the distance scaled
# by 20000/100000. No shipped map has a gap in that range.
awk '{ if (NR > 4162) { $4 += 100000; $3 = sprintf("%.15g", $3 + 0.1) } print }' example.map > gap100.map

# fd1/fd2.vcf: p1/p2 with every 400th site (from the 200th) turned into a fixed
# difference, all 0|0 in p1 and all 1|1 in p2. With equal panel sizes these
# sites sit at joint frequency exactly 0.5 and carry no heterozygotes.
for k in 1 2; do
    gt="0|0"; [ $k = 2 ] && gt="1|1"
    gzip -dc p$k.vcf.gz | awk -v gt="$gt" 'BEGIN{OFS="\t"} /^#/ {print; next}
        { r++; if (r % 400 == 200) { for (i=10;i<=NF;i++) $i=gt; print $3 > "fd_ids.txt" } print }' > fd$k.vcf
done

# c1/c2.tped: p1/p2 converted from VCF to TPED, record for record, to check
# that TPED input gives exactly what the same data gives as VCF. Nothing is
# removed: duplicated positions (12 in p1) and, in XP mode, sites monomorphic
# in the two panels combined are dropped by every loader.
for k in 1 2; do
    gzip -dc p$k.vcf.gz > dd$k.vcf
    awk '/^#/ {next} { printf "%s %s %.10g %s", $1, $3, $2/1e6, $2;
        for (i=10;i<=NF;i++) printf " %s %s", substr($i,1,1), substr($i,3,1); printf "\n" }' dd$k.vcf > c$k.tped
done

# dup.map: example.map with Locus2000 moved onto Locus1999's position (map rows
# 2000-2001, both reported as iHS cores on the normal map), so two records
# share a position; both must be dropped from hap/thap input.
awk 'NR==2000 {p=$4} { if (NR==2001) $4=p; print }' example.map > dup.map

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
golden ehh12-alt      ehh12.Locus4077.alt.out -- --ehh12 Locus4077 --alt --hap example.hap --map example.map
golden ehh12-minor    ehh12.Locus1003.out -- --ehh12 Locus1003 --hap example.hap --map example.map
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
golden ihs-gap100     ihs.out          -- --ihs --hap example.hap --map gap100.map
golden nsl-gap100     nsl.out          -- --nsl --hap example.hap --map gap100.map
golden ihh12-gap100   ihh12.out        -- --ihh12 --hap example.hap --map gap100.map

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
    # Statistics that use every site (pi, EHH, XP) combined with ones that
    # MAF-filter (iHS): each side must still equal its single run.
    run s.log --ihs --vcf p1.vcf.gz --pmap --out s_ihs
    run s.log --pi --vcf p1.vcf.gz --pmap --out s_pi
    run s.log --ehh chr1:29243 --vcf p1.vcf.gz --pmap --out s_ehh
    run s.log --xpehh --vcf p1.vcf.gz --vcf-ref p2.vcf.gz --pmap --out s_xp
    run m.log --ihs --pi --vcf p1.vcf.gz --pmap --out m_pi
    identical "multistat==single:pi-with-ihs" m_pi.pi.100bp.out s_pi.pi.100bp.out
    identical "multistat==single:ihs-with-pi" m_pi.ihs.out s_ihs.ihs.out
    run m.log --ihs --ehh chr1:29243 --vcf p1.vcf.gz --pmap --out m_ehh
    identical "multistat==single:ehh-with-ihs" m_ehh.ehh.chr1:29243.out s_ehh.ehh.chr1:29243.out
    run m.log --ihs --xpehh --vcf p1.vcf.gz --vcf-ref p2.vcf.gz --pmap --out m_xp
    identical "multistat==single:xpehh-with-ihs" m_xp.xpehh.out s_xp.xpehh.out
    identical "multistat==single:ihs-with-xpehh" m_xp.ihs.out s_ihs.ihs.out
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
    # ... and the same records as thap + map (one line per locus)
    for k in 1 2; do
        awk '/^#/ {next} { s = ""; for (i=10;i<=NF;i++) s = s " " substr($i,1,1) " " substr($i,3,1); print substr(s,2) }' dd$k.vcf > c$k.thap
    done
    awk '/^#/ {next} { printf "%s %s %.10g %s\n", $1, $3, $2/1e6, $2 }' dd1.vcf > c.map
    run te5.log --xpehh --thap c1.thap --thap-ref c2.thap --map c.map --pmap --out te_h
    identical "thap-equiv:xpehh" te_h.xpehh.out te_v.xpehh.out
    # hap input with a duplicated map position: both records are dropped and logged
    run dp.log --ihs --hap example.hap --map dup.map --out dp
    nd=$(awk -F'\t' 'NR>1 && ($2=="Locus1999" || $2=="Locus2000")' dp.ihs.out | wc -l | tr -d ' ')
    if [ "$RC" -eq 0 ] && [ "$nd" = 0 ] && grep -q "Dropped 2 records at 1 duplicated positions" dp.log; then pass "dup-positions:hap"
    else fail "dup-positions:hap" "exit $RC, $nd of the duplicated loci reported, log note $(grep -q 'Dropped 2 records' dp.log && echo present || echo absent)"; fi
fi

if selected xp-intersect; then
    # XP VCF input is paired by position and REF/ALT (B15). Start from p1/p2
    # without duplicated positions (u1/u2 list the same positions), then
    # remove every 13th record from file 1, every 7th from file 2, and change
    # ALT on every 11th record of file 2. The result must equal a run on a
    # matched pair with the union of those records removed from both files.
    for k in 1 2; do
        awk 'NR==FNR {if (!/^#/) c[$2]++; next} /^#/ || c[$2]==1' dd$k.vcf dd$k.vcf > xu$k.vcf
        awk '/^#/ {print; next} {n++} n%13 && n%7 && n%11' xu$k.vcf > xr$k.vcf
    done
    awk '/^#/ {print; next} {n++} n%13 != 0' xu1.vcf > xq1.vcf
    awk -F'\t' -v OFS='\t' '/^#/ {print; next} {n++} n%7 == 0 {next} n%11 == 0 {$5 = ($5=="A" ? "C" : "A")} {print}' xu2.vcf > xq2.vcf
    run xq.log --xpehh --vcf xq1.vcf --vcf-ref xq2.vcf --pmap --out xq
    run xr.log --xpehh --vcf xr1.vcf --vcf-ref xr2.vcf --pmap --out xr
    identical "xp-intersect:vcf" xq.xpehh.out xr.xpehh.out
    # XP-nSL counts retained sites only (B47), so records removed before the
    # analysis must not change it either
    for ph in "" "--unphased"; do
        run xqn.log --xpnsl $ph --vcf xq1.vcf --vcf-ref xq2.vcf --pmap --out xqn
        run xrn.log --xpnsl $ph --vcf xr1.vcf --vcf-ref xr2.vcf --pmap --out xrn
        identical "xp-intersect:xpnsl${ph:+-unphased}" xqn.xpnsl.out xrn.xpnsl.out
    done
    # expected drop counts, computed from the record numbers
    nrec=$(grep -vc '^#' xu1.vcf)
    exp=$(awk -v N="$nrec" 'BEGIN {for (n=1;n<=N;n++) {a+=(n%7==0 && n%13); b+=(n%13==0 && n%7); m+=(n%11==0 && n%7 && n%13)} print a, b, m}')
    got="$(sed -n 's/^Dropped \([0-9]*\) records of xq1.vcf at positions absent.*/\1/p' xq.log) $(sed -n 's/^Dropped \([0-9]*\) records of xq2.vcf at positions absent.*/\1/p' xq.log) $(sed -n 's/^Dropped \([0-9]*\) shared positions where REF\/ALT differ.*/\1/p' xq.log)"
    if [ "$got" = "$exp" ]; then pass "xp-intersect:log-counts"
    else fail "xp-intersect:log-counts" "logged '$got', expected '$exp'"; fi
    # Equal record counts but offset by one record: 3.0.1 paired records by
    # index and silently reported every site against the wrong partner.
    awk -v N="$nrec" '/^#/ {print; next} {n++} n<N' xu1.vcf > xs1.vcf
    awk '/^#/ {print; next} {n++} n>1' xu2.vcf > xs2.vcf
    for k in 1 2; do awk -v N="$nrec" '/^#/ {print; next} {n++} n>1 && n<N' xu$k.vcf > xm$k.vcf; done
    run xs.log --xpehh --vcf xs1.vcf --vcf-ref xs2.vcf --pmap --out xs
    run xm.log --xpehh --vcf xm1.vcf --vcf-ref xm2.vcf --pmap --out xm
    identical "xp-intersect:offset" xs.xpehh.out xm.xpehh.out
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

if selected maf-symmetry; then
    # The MAF filter is symmetric: swapping 0 and 1 at every site (derived
    # frequency f -> 1-f) must keep exactly the same cores. At --maf 0.1 with
    # 20 haplotypes, sites with 2 and with 18 derived alleles are both exactly
    # at the boundary and both kept (ties are kept).
    gzip -dc p1.vcf.gz | awk 'BEGIN { OFS = "\t" } /^#/ { print; next }
        { for (i = 10; i <= NF; i++) { a = substr($i,1,1); b = substr($i,3,1)
              $i = (a == "0" ? "1" : "0") substr($i,2,1) (b == "0" ? "1" : "0") } print }' > flip.vcf
    run ms1.log --ihs --vcf p1.vcf.gz --pmap --maf 0.1 --out ms_a
    run ms2.log --ihs --vcf flip.vcf --pmap --maf 0.1 --out ms_b
    r=$(awk -F'\t' 'FNR == 1 { next } FNR == NR { a[$2] = 1; na++; next } { nb++; if (!($2 in a)) only++ }
        END { if (na != nb || only) printf "%d cores vs %d after flipping alleles (%d not shared)", na, nb, only }' ms_a.ihs.out ms_b.ihs.out)
    if [ -s ms_a.ihs.out ] && [ -z "$r" ]; then pass "maf-symmetry:ihs"; else fail "maf-symmetry:ihs" "${r:-no output}"; fi
fi

if selected pi-edges; then
    # pi-large-n: 23,200 diploids = 46,400 haplotypes, above the 46,342 at which
    # nhaps*(nhaps-1) overflowed int. Site 2 carries the derived allele on
    # exactly half the haplotypes, so window 101-200 has pi = (n/2)^2/(n(n-1)/2).
    awk 'BEGIN { n = 23200
        print "##fileformat=VCFv4.2"
        print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
        for (i = 1; i <= n; i++) printf "\ts%d", i
        printf "\n"
        for (s = 1; s <= 3; s++) {
            printf "1\t%d\tsnp%d\tA\tC\t.\tPASS\t.\tGT", s * 75, s
            for (i = 1; i <= n; i++) printf "\t%s", (s == 2 ? "0|1" : (i <= n / 4 ? "1|1" : "0|0"))
            printf "\n"
        } }' > pi_large.vcf
    run pl.log --pi --pi-win 100 --vcf pi_large.vcf --pmap --out pl
    r=$(awk -v n=46400 '
        NR == 2 { want = (n/2) * (n/2) / (n * (n-1) / 2); d = $3 - want; if (d < 0) d = -d
                  if (d > 1e-5 * want) printf "window 101-200 pi = %s, expected %.6g; ", $3, want }
        END { if (NR != 3) printf "%d windows, expected 3", NR }' pl.pi.100bp.out 2>/dev/null)
    if [ "$RC" -eq 0 ] && [ -s pl.pi.100bp.out ] && [ -z "$r" ]; then pass "pi-large-n"
    else fail "pi-large-n" "${r:-exit $RC}"; fi
    # pi-edge: the last position is 0 and then -1 (mod --pi-win); both read past
    # the last locus. The window list must end with the window holding it.
    for last in 300 299; do
        gzip -dc p1.vcf.gz | awk -v last=$last 'BEGIN { OFS = "\t" } /^#/ { print; next }
            { r++; if (r <= 3) { $2 = (r == 3 ? last : r * 50); print } }' > pe.vcf
        run pe.log --pi --pi-win 100 --vcf pe.vcf --pmap --out pe
        nw=0; [ -f pe.pi.100bp.out ] && nw=$(($(wc -l < pe.pi.100bp.out)))
        if [ "$RC" -eq 0 ] && [ "$nw" = 3 ]; then pass "pi-edge:last=$last"
        else fail "pi-edge:last=$last" "exit $RC, $nw windows (expected 3)"; fi
        rm -f pe.pi.100bp.out
    done
fi

if selected xp-validate; then
    # B16: the XP VCF reader validates genotypes like the single-population
    # reader. Record 104 of p2 (position 16202, first genotype 1|1) is
    # corrupted (missing allele, '/' in phased mode, a genotype field
    # removed) or the #CHROM line is deleted; XP must
    # exit 1, as --ihs does on the same file. 3.0.1 (7a93ff1) runs the first
    # three with exit 0; on the missing header it prints "number of
    # haplotypes (-2) and number of loci (1483) must be positive" and aborts
    # (exit 134).
    awk -F'\t' -v OFS='\t' '/^#/ {print; next} {n++} n==104 {$10=".|."} {print}' dd2.vcf > xv_miss.vcf
    awk -F'\t' -v OFS='\t' '/^#/ {print; next} {n++} n==104 {$10="0/1"} {print}' dd2.vcf > xv_slash.vcf
    awk -F'\t' -v OFS='\t' '/^#/ {print; next} {n++} n==104 {NF=NF-1} {print}' dd2.vcf > xv_short.vcf
    grep -v '^#CHROM' dd2.vcf > xv_nohdr.vcf
    for c in miss slash short nohdr; do
        run xv.log --xpehh --vcf dd1.vcf --vcf-ref xv_$c.vcf --pmap --out xv; xrc=$RC
        run xvi.log --ihs --vcf xv_$c.vcf --pmap --out xvi
        if [ "$xrc" -eq 1 ] && [ "$RC" -eq 1 ]; then pass "xp-validate:$c"
        else fail "xp-validate:$c" "xpehh exit $xrc, ihs exit $RC (expected 1 and 1)"; fi
    done
    # With --filter-vcf a bad record is skipped instead, and its partner in the
    # other file goes with it. The record is polymorphic, so reading '.' as 0
    # (3.0.1) changes the result (REF/ALT relabelled A/C so --filter-vcf keeps SNPs)
    ac='/^#/ {print; next} {n++; $4="A"; $5="C"} '
    awk -F'\t' -v OFS='\t' "$ac"'{print}' dd1.vcf > xv_a1.vcf
    awk -F'\t' -v OFS='\t' "$ac"'{print}' xv_miss.vcf > xv_a2.vcf
    for k in 1 2; do awk -F'\t' -v OFS='\t' "$ac"'n != 104 {print}' dd$k.vcf > xv_r$k.vcf; done
    run xva.log --xpehh --filter-vcf --vcf xv_a1.vcf --vcf-ref xv_a2.vcf --pmap --out xva
    run xvr.log --xpehh --filter-vcf --vcf xv_r1.vcf --vcf-ref xv_r2.vcf --pmap --out xvr
    identical "xp-validate:filter-vcf-skip" xva.xpehh.out xvr.xpehh.out
fi

if selected norm-robust; then
    # Inputs: two iHS tables from the goldens (ihs: p1 VCF; ihs-hap:
    # example.hap) and six 60 kb genes across the p1 range (49-481 kb).
    gzip -dc "$EXPECTED/ihs.tsv.gz" > nr_a.ihs.out
    gzip -dc "$EXPECTED/ihs-hap.tsv.gz" > nr_b.ihs.out
    awk 'BEGIN { for (g=0; g<6; g++) printf "chr1\t%d\t%d\tG%d\n", 50000+g*70000, 110000+g*70000, g+1 }' > nr_genes.bed
    run nr0.log norm --ihs --bins 20 --files nr_a.ihs.out
    # B23: --norm-files + --gene-bed without --bp-win deleted uninitialised
    # pointers after writing its output (3.0.1: exit 134/139)
    run nr1.log norm --ihs --norm-files nr_a.ihs.out.20bins.norm --gene-bed nr_genes.bed --min-snps 1
    if [ "$RC" -eq 0 ] && [ -s nr_genes.bed.ihs.genetable ]; then pass "norm-gene:norm-files"
    else fail "norm-gene:norm-files" "exit $RC"; fi
    rm -f nr_genes.bed.ihs.genetable
    # B24: --files a b --first --gene-bed wrote past the end of a vector. The
    # release build hides it (undefined behaviour); a libc++-hardened build
    # aborted. Regression check only.
    run nr2.log norm --ihs --bins 20 --files nr_a.ihs.out nr_b.ihs.out --first --gene-bed nr_genes.bed --min-snps 1
    if [ "$RC" -eq 0 ] && [ -s nr_genes.bed.ihs.genetable ]; then pass "norm-gene:first"
    else fail "norm-gene:first" "exit $RC"; fi
    rm -f nr_genes.bed.ihs.genetable
    # B25: with fewer than 3 genes the length regression is skipped and the
    # adjusted score is NA; 3.0.1 aborted in GSL with one gene
    head -1 nr_genes.bed > nr_one.bed
    run nr3.log norm --ihs --norm-files nr_a.ihs.out.20bins.norm --gene-bed nr_one.bed --min-snps 1
    adj=$(awk -F'\t' 'NR==2 {print $9}' nr_one.bed.ihs.genetable 2>/dev/null)
    if [ "$RC" -eq 0 ] && [ "$adj" = NA ] && grep -q 'length regression' nr3.log; then pass "norm-gene:one-gene"
    else fail "norm-gene:one-gene" "exit $RC, top_score_adj '$adj'"; fi
    # B26: XP bottom_score_adj must be the residual of bottom_score on
    # log(len), recomputed here by least squares (genes of unequal length).
    gzip -dc "$EXPECTED/xpehh.tsv.gz" > nr_x.xpehh.out
    run nr4.log norm --xpehh --files nr_x.xpehh.out
    awk 'BEGIN { s=50000; for (g=0; g<6; g++) { e=s+30000+g*12000; printf "chr1\t%d\t%d\tX%d\n", s, e, g+1; s=e+10000 } }' > nr_xg.bed
    run nr5.log norm --xpehh --norm-files nr_x.xpehh.out.norm --gene-bed nr_xg.bed --min-snps 1
    d=$(awk -F'\t' 'NR>1 { n++; L[n]=log($5); y[n]=$11; a[n]=$12; sl+=L[n]; sy+=y[n] }
        END { if (n < 3) { print "n=" n; exit } ml=sl/n; my=sy/n
              for (i=1;i<=n;i++) { sxy+=(L[i]-ml)*(y[i]-my); sxx+=(L[i]-ml)^2 }
              b=sxy/sxx; m=0; for (i=1;i<=n;i++) { r=y[i]-my-b*(L[i]-ml)-a[i]; if (r<0) r=-r; if (r>m) m=r }
              printf "%.3g", m }' nr_xg.bed.xpehh.genetable 2>/dev/null)
    if [ "$RC" -eq 0 ] && awk -v d="$d" 'BEGIN { exit !(d != "" && d+0 < 1e-4) }'; then pass "norm-gene:xp-bottom-adj"
    else fail "norm-gene:xp-bottom-adj" "exit $RC, max |bottom_score_adj - residual| = '$d'"; fi
    # B37: BED genes cover 1-based start+1..end inclusive. A gene ending
    # exactly at a SNP must count it (3.0.1 dropped the last base).
    P=$(awk -F'\t' 'NR==202 {print $3}' nr_a.ihs.out.20bins.norm)
    nb=$(awk -F'\t' -v P="$P" 'NR>1 && $3 >= P-999 && $3 <= P' nr_a.ihs.out.20bins.norm | wc -l | tr -d ' ')
    printf 'chr1\t%d\t%d\tENDP\nchr1\t100000\t150000\tF1\nchr1\t300000\t350000\tF2\n' $((P-1000)) "$P" > nr_end.bed
    run nr6.log norm --ihs --norm-files nr_a.ihs.out.20bins.norm --gene-bed nr_end.bed --min-snps 1
    ns=$(awk -F'\t' '$4=="ENDP" {print $6}' nr_end.bed.ihs.genetable 2>/dev/null)
    if [ "$RC" -eq 0 ] && [ "$ns" = "$nb" ]; then pass "norm-gene:snp-at-gene-end"
    else fail "norm-gene:snp-at-gene-end" "exit $RC, nsnps '$ns', expected $nb"; fi
    # ... and windows (1-based inclusive) overlap a gene that starts on their
    # last base or ends on their first base (--win-files path)
    printf 'chr\tstart\tend\tn_snps\tfrac_extreme\tperc\tscore\nchr1\t1\t100000\t5\t0\t1\t1\nchr1\t100001\t200000\t5\t0\t1\t1\nchr1\t200001\t300000\t5\t0\t1\t1\n' > nr_w.windows
    printf 'chr1\t199999\t260000\tSTARTS\nchr1\t50000\t100001\tENDS\n' > nr_wg.bed
    run nr7.log norm --ihs --win-files nr_w.windows --gene-bed nr_wg.bed
    ov=$(awk -F'\t' '$2==100001 {print $8}' nr_w.windows.ann 2>/dev/null)
    if [ "$RC" -eq 0 ] && [ "$ov" = "ENDS, STARTS" ]; then pass "norm-gene:window-edges"
    else fail "norm-gene:window-edges" "exit $RC, window 100001-200000 overlaps '$ov' (expected 'ENDS, STARTS')"; fi
    # B48: --bp-win --gene-bed must annotate the .windows table it wrote, as
    # --win-files does; 3.0.1 read the .norm file as windows and wrote
    # .norm.ann with every SNP tagged
    cp nr_a.ihs.out.20bins.norm nr_bw.norm
    run nr8.log norm --ihs --norm-files nr_bw.norm --bp-win --winsize 100000 --gene-bed nr_wg.bed --min-snps 1
    mv nr_bw.norm.100kb.windows.ann nr_bw_direct.ann 2>/dev/null
    run nr9.log norm --ihs --win-files nr_bw.norm.100kb.windows --gene-bed nr_wg.bed
    identical "norm-gene:bpwin-annotates-windows" nr_bw_direct.ann nr_bw.norm.100kb.windows.ann
    # B38: the permutation test reads norm's own gene tables (top_score_adj),
    # is reproducible with --perm-seed, and reports (count+1)/(n+1). The
    # six-gene XP table is split into alternate genes.
    awk 'NR==1 || NR%2==0' nr_xg.bed.xpehh.genetable > nr_setA.tab
    awk 'NR==1 || NR%2==1' nr_xg.bed.xpehh.genetable > nr_setB.tab
    run nr10.log norm --gene-target nr_setA.tab --gene-background nr_setB.tab --perm-seed 11; r1=$RC
    p1=$(sed -n 's/.*one-sided p-value (A > B): //p' nr10.log)
    run nr11.log norm --gene-target nr_setA.tab --gene-background nr_setB.tab --perm-seed 11
    p2=$(sed -n 's/.*one-sided p-value (A > B): //p' nr11.log)
    if [ "$r1" -eq 0 ] && [ "$RC" -eq 0 ] && [ -n "$p1" ] && [ "$p1" = "$p2" ] && awk -v p="$p1" 'BEGIN { exit !(p > 0 && p <= 1) }'
    then pass "norm-gene:perm-test"
    else fail "norm-gene:perm-test" "exit $r1/$RC, p-values '$p1' '$p2'"; fi
    # B18: a non-numeric or non-finite score is treated as missing: the output
    # equals a run with that row deleted, and the count is logged. 3.0.1's
    # stream entered a fail state (iHS: 33 of 527 rows written, exit 0).
    for c in "ihs|ihs|7|nan|20bins.norm|--ihs --bins 20" "xpehh|xpehh|9|x|norm|--xpehh" "ihh12|ihh12|5|inf|norm|--ihh12"; do
        IFS='|' read -r st g col tok suf args <<< "$c"
        gzip -dc "$EXPECTED/$g.tsv.gz" | awk -F'\t' -v OFS='\t' -v c="$col" -v t="$tok" 'NR==101 {$c=t} {print}' > nb_bad.$st.out
        gzip -dc "$EXPECTED/$g.tsv.gz" | awk 'NR!=101' > nb_rm.$st.out
        run nb1.log norm $args --files nb_bad.$st.out
        grep -q "WARNING: 1 rows have a non-numeric" nb1.log; w=$?
        run nb2.log norm $args --files nb_rm.$st.out
        if [ "$w" -ne 0 ]; then fail "norm-bad-token:$st" "no warning logged for '$tok'"
        else identical "norm-bad-token:$st" nb_bad.$st.out.$suf nb_rm.$st.out.$suf; fi
    done
    # B19: --crit-percent cutoffs come from finite normalised scores only, so
    # adding a missing (-9999) locus must not move them (3.0.1 sorted -9999
    # into the lower tail).
    gzip -dc "$EXPECTED/ihs-hap.tsv.gz" > nc_a.ihs.out
    cp nc_a.ihs.out nc_b.ihs.out
    printf 'chr1\tmissing1\t9999999\t0.5\t0.1\t0.1\t-9999\n' >> nc_b.ihs.out
    run nc1.log norm --ihs --crit-percent 0.02 --files nc_a.ihs.out
    c1=$(grep -E '^(Top|Bottom) cutoff' nc1.log | tr '\n' ' ')
    run nc2.log norm --ihs --crit-percent 0.02 --files nc_b.ihs.out
    c2=$(grep -E '^(Top|Bottom) cutoff' nc2.log | tr '\n' ' ')
    if [ -n "$c1" ] && [ "$c1" = "$c2" ]; then pass "norm-critpct:missing-ignored"
    else fail "norm-critpct:missing-ignored" "cutoffs '$c1' vs '$c2'"; fi
    # U10: iHH12 is one-tailed; --crit-percent 0.05 flags about 5% of loci
    # (3.0.1 used 0.025 per tail and only the top tail is flagged)
    gzip -dc "$EXPECTED/ihh12.tsv.gz" > nc.ihh12.out
    run nc3.log norm --ihh12 --crit-percent 0.05 --files nc.ihh12.out
    fr=$(awk -F'\t' 'NR>1 {n++; f+=($NF==1)} END {if (n) printf "%.4f", f/n}' nc.ihh12.out.norm 2>/dev/null)
    if [ "$RC" -eq 0 ] && awk -v f="$fr" 'BEGIN { exit !(f != "" && f >= 0.045 && f <= 0.055) }'; then pass "norm-critpct:ihh12-one-tailed"
    else fail "norm-critpct:ihh12-one-tailed" "exit $RC, fraction flagged '$fr' (expected 0.05)"; fi
    # B21: --bp-win writes the window holding each file's last SNP (3.0.1
    # never closed it). The last window row must contain the last position.
    for c in "ihs|ihs-hap|--ihs --bins 20|20bins.norm" "xpehh|xpehh|--xpehh|norm" "ihh12|ihh12|--ihh12|norm"; do
        IFS='|' read -r st g args suf <<< "$c"
        gzip -dc "$EXPECTED/$g.tsv.gz" > nw.$st.out
        run nw.log norm $args --bp-win --winsize 100000 --files nw.$st.out
        last=$(awk -F'\t' 'NR>1 && $3>m {m=$3} END {print m}' nw.$st.out.$suf 2>/dev/null)
        wl=$(tail -1 nw.$st.out.$suf.100kb.windows 2>/dev/null | awk -F'\t' '{print $2"-"$3}')
        if [ "$RC" -eq 0 ] && [ -n "$last" ] && awk -v w="$wl" -v p="$last" 'BEGIN { split(w,a,"-"); exit !(p+0 >= a[1]+0 && p+0 <= a[2]+0) }'
        then pass "norm-bpwin:last-window:$st"
        else fail "norm-bpwin:last-window:$st" "exit $RC, last SNP $last, last window $wl"; fi
    done
    # B45: iHH12 window maxima are per file. File A ends in a large score;
    # file B's first window must report B's own maximum.
    gzip -dc "$EXPECTED/ihh12.tsv.gz" > nm_B.ihh12.out
    awk -F'\t' -v OFS='\t' -v n="$(wc -l < nm_B.ihh12.out)" 'NR==n {$5=50} {print}' nm_B.ihh12.out > nm_A.ihh12.out
    run nm.log norm --ihh12 --bp-win --winsize 100000 --files nm_A.ihh12.out nm_B.ihh12.out
    w1=$(awk -F'\t' 'NR==2 {print $2, $3, $7}' nm_B.ihh12.out.norm.100kb.windows 2>/dev/null)
    mx=$(echo "$w1" | { read s e sc; awk -F'\t' -v s="$s" -v e="$e" 'NR>1 && $3>=s && $3<=e && (m=="" || $6>m) {m=$6} END {print m}' nm_B.ihh12.out.norm; })
    sc=$(echo "$w1" | awk '{print $3}')
    if [ "$RC" -eq 0 ] && [ -n "$mx" ] && awk -v a="$sc" -v b="$mx" 'BEGIN { d=a-b; if (d<0) d=-d; exit !(d < 1e-4) }'
    then pass "norm-bpwin:ihh12-per-file-max"
    else fail "norm-bpwin:ihh12-per-file-max" "exit $RC, first window of file B scored '$sc', its SNPs' maximum '$mx'"; fi
    # B22/U4: 995 windows with 10 SNPs (every 20th has one extreme SNP) and 5
    # with 12 SNPs and none. Tied SNP-count bounds left a bin unfilled in
    # 3.0.1, whose windows were compared with 0 and labelled 1%; with '>='
    # windows tied at a 0 boundary were labelled too. No window without an
    # extreme SNP may be labelled 1 or 5.
    awk 'BEGIN { OFS="\t"; print "chr","id","pos","freq","ihh1","ihh0","ihs","norm_ihs","crit"
        for (k=0; k<1000; k++) { n = (k < 995) ? 10 : 12
            for (j=0; j<n; j++) { p = 1 + k*100000 + j*1000; c = (k < 995 && k%20==0 && j==0)
                print "chr1", "chr1:" p, p, 0.5, 0.1, 0.1, 0, (c ? 3 : 0.1), c } } }' > nq.ihs.out.20bins.norm
    run nq.log norm --ihs --norm-files nq.ihs.out.20bins.norm --bp-win --winsize 100000
    bad=$(awk -F'\t' 'NR>1 && $5==0 && ($6==1 || $6==5)' nq.ihs.out.20bins.norm.100kb.windows 2>/dev/null | wc -l | tr -d ' ')
    nw=$(awk 'END {print NR-1}' nq.ihs.out.20bins.norm.100kb.windows 2>/dev/null)
    if [ "$RC" -eq 0 ] && [ "$nw" = 1000 ] && [ "$bad" = 0 ]; then pass "norm-bpwin:quantile-ties"
    else fail "norm-bpwin:quantile-ties" "exit $RC, $nw windows, $bad windows without extreme SNPs labelled 1 or 5"; fi
    # B49: the --log file records the window analysis in both modes (3.0.1:
    # closed before it with --files, never written with --norm-files)
    rm -f nl1.log nl2.log
    run nlr1.txt norm --ihs --bins 20 --bp-win --winsize 100000 --files nr_b.ihs.out --log nl1.log
    run nlr2.txt norm --ihs --norm-files nr_b.ihs.out.20bins.norm --bp-win --winsize 100000 --log nl2.log
    for k in 1 2; do
        if grep -q 'Creating window file' nl$k.log 2>/dev/null; then pass "norm-log:window-analysis:$k"
        else fail "norm-log:window-analysis:$k" "window analysis missing from the --log file"; fi
    done
    # B20 (Q11): --log-input takes the bin count from the log and reproduces
    # the run that wrote it exactly (statistics are logged at full
    # precision). 3.0.1 always used 100 bins and silently dropped loci.
    gzip -dc "$EXPECTED/ihs.tsv.gz" > li_ref.ihs.out; cp li_ref.ihs.out li_use.ihs.out
    run li1.log norm --ihs --bins 20 --files li_ref.ihs.out --log li_ref.log
    run li2.log norm --ihs --files li_use.ihs.out --log-input li_ref.log --log li_use.log
    identical "norm-log-input:ihs-20-bins" li_use.ihs.out.20bins.norm li_ref.ihs.out.20bins.norm
    gzip -dc "$EXPECTED/xpehh.tsv.gz" > li_ref.xpehh.out; cp li_ref.xpehh.out li_use.xpehh.out
    run li3.log norm --xpehh --files li_ref.xpehh.out --log li_refx.log
    run li4.log norm --xpehh --files li_use.xpehh.out --log-input li_refx.log --log li_usex.log
    identical "norm-log-input:xpehh" li_use.xpehh.out.norm li_ref.xpehh.out.norm
    run li5.log norm --ihs --bins 30 --files li_use.ihs.out --log-input li_ref.log --log li_bad.log
    if [ "$RC" -eq 1 ] && grep -q 'disagrees' li5.log; then pass "norm-log-input:bins-conflict"
    else fail "norm-log-input:bins-conflict" "exit $RC (expected 1 with a --bins conflict message)"; fi
    # U2: --log-only writes the bin statistics and no normalised output, for
    # every statistic; window and gene analyses are skipped (3.0.1 ignored
    # it for iHS/nSL and aborted on XP with --bp-win)
    gzip -dc "$EXPECTED/ihs.tsv.gz" > lo.ihs.out; gzip -dc "$EXPECTED/xpehh.tsv.gz" > lo.xpehh.out
    run lo1.log norm --ihs --log-only --bp-win --files lo.ihs.out --log lo1.norm.log; r1=$RC
    run lo2.log norm --xpehh --log-only --bp-win --files lo.xpehh.out --log lo2.norm.log; r2=$RC
    nout=$(ls lo.ihs.out.*norm* lo.xpehh.out.*norm* 2>/dev/null | wc -l | tr -d ' ')
    if [ "$r1" -eq 0 ] && [ "$r2" -eq 0 ] && [ "$nout" = 0 ] && grep -q 'mean' lo1.norm.log; then pass "norm-log-only"
    else fail "norm-log-only" "exit $r1/$r2, $nout normalised or window files written"; fi
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
