#!/usr/bin/env bash
# =============================================================================
# verify_and_stage_beds.sh - smoke-test and finalise the ClinCNV BED set
# =============================================================================
# Usage:
#   ./verify_and_stage_beds.sh --bam /path/to/a/sample.bam [--only v5,v6]
#
# Runs BedCoverage against each staged BED with a real BAM, checks the output
# column count, then writes checksums and a provenance note.
# =============================================================================

set -uo pipefail

BEDDIR="${BEDDIR:-path/to/reference/EXOME_XtraXXX/clincnv_beds}"
SIF_DIR="${SIF_DIR:-path/to/Docker_Images}"
BUILD_SCRIPT="${BUILD_SCRIPT:-path/to/reference/EXOME_XtraXXX/build_native_bed.sh}"
BAM=""
ONLY=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam)    BAM="$2"; shift 2 ;;
    --beddir) BEDDIR="$2"; shift 2 ;;
    --only)   ONLY="$2"; shift 2 ;;
    -h|--help) sed -n '2,12p' "$0"; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

NGSBITS="$SIF_DIR/ngsbits.sif"
hr() { printf '%s\n' "----------------------------------------------------------"; }

cd "$BEDDIR" || { echo "ERROR: cannot cd to $BEDDIR" >&2; exit 1; }

# select which BEDs to check (all by default; --only v5,v6 to subset)
BEDS=()
for f in *.bed; do
  [[ -e "$f" ]] || continue
  if [[ -n "$ONLY" ]]; then
    keep=0
    IFS=',' read -ra pats <<< "$ONLY"
    for pat in "${pats[@]}"; do [[ "$f" == *"$pat"* ]] && keep=1; done
    [[ "$keep" -eq 1 ]] || continue
  fi
  BEDS+=("$f")
done
[[ ${#BEDS[@]} -gt 0 ]] || { echo "ERROR: no BEDs matched in $BEDDIR" >&2; exit 1; }

echo "=== BED SET ==="
[[ -n "$ONLY" ]] && echo "(filtered: $ONLY)"
printf '%s\n' "${BEDS[@]}"
hr

# --- structural checks on every BED ----------------------------------------
echo "=== STRUCTURE ==="
printf "%-42s %9s %13s %6s %s\n" FILE REGIONS BASES COLS CHRORDER
for f in "${BEDS[@]}"; do
  cols=$(awk '{print NF}' "$f" | sort -u | tr '\n' ',' | sed 's/,$//')
  regions=$(wc -l < "$f")
  bases=$(awk '{s+=$3-$2} END{print s+0}' "$f")
  # ngs-bits order: chr1..chr22, chrX, chrY
  order=$(cut -f1 "$f" | uniq | tr '\n' ' ')
  expect="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY "
  [[ "$order" == "$expect" ]] && ok="OK" || ok="MISMATCH"
  printf "%-42s %9s %13s %6s %s\n" "$f" "$regions" "$bases" "$cols" "$ok"
  [[ "$ok" == "MISMATCH" ]] && echo "    got: $order"
  # GC sanity
  bad=$(awk '$4 !~ /^[0-9]*\.?[0-9]+$/ || $4<0 || $4>1' "$f" | wc -l)
  [[ "$bad" -gt 0 ]] && echo "    WARNING: $bad rows with invalid GC in col4"
done
hr

# --- BedCoverage smoke test -------------------------------------------------
echo "=== BEDCOVERAGE SMOKE TEST ==="
if [[ -z "$BAM" ]]; then
  echo "SKIPPED: no --bam supplied."
elif [[ ! -f "$BAM" ]]; then
  echo "ERROR: BAM not found: $BAM"
else
  echo "bam: $BAM"
  BAMDIR=$(dirname "$BAM")
  for f in "${BEDS[@]}"; do
    printf "  %-42s " "$f"
    out=$(singularity exec \
            -B /home/data/human_genetics -B /home/sdata -B "$BAMDIR" \
            "$NGSBITS" BedCoverage -bam "$BAM" -in "$BEDDIR/$f" \
            -min_mapq 5 -decimals 4 2>"$BEDDIR/.bedcov.err" | head -3)
    if [[ -z "$out" ]]; then
      echo "FAILED"
      sed 's/^/      /' "$BEDDIR/.bedcov.err" | head -5
    else
      ncol=$(echo "$out" | sed -n '2p' | awk '{print NF}')
      if [[ "$ncol" -eq 6 ]]; then
        echo "OK (${ncol} cols -> strip f4,f5 leaves chr/start/end/coverage)"
      else
        echo "UNEXPECTED: ${ncol} cols"
        echo "$out" | sed 's/^/      /'
      fi
    fi
  done
  rm -f "$BEDDIR/.bedcov.err"
fi
hr

# --- stage build script and checksums --------------------------------------
echo "=== STAGING ==="
if [[ -f "$BUILD_SCRIPT" ]] && [[ ! -f "$BEDDIR/$(basename "$BUILD_SCRIPT")" ]]; then
  cp "$BUILD_SCRIPT" "$BEDDIR/" && echo "copied $(basename "$BUILD_SCRIPT")"
fi

md5sum *.bed > checksums.txt   # all BEDs, not just the filtered set
echo "wrote checksums.txt:"
sed 's/^/  /' checksums.txt

# --- provenance note --------------------------------------------------------
if [[ -z "$ONLY" || -f PROVENANCE.md ]]; then :; fi
cat > PROVENANCE.md << 'EOF'
# ClinCNV target designs

Built for germline CNV calling. One BED per capture-kit group; samples are
never pooled across groups, because ClinCNV normalises read depth against the
cohort and a region absent from a sample's capture design reads as zero depth
(i.e. a false homozygous deletion).

## Files

| BED | Group | Source designs | Samples |
|---|---|---|---|
| legacy_consensus.min20.annotated.bed | legacy | S000032/S000042 (v1/v2) + S000048/S000052 (v3/v4) | 71 |
| v5_native.min20.annotated.bed | v5 | S000054 cexome5 CNV targets | 39 |
| v6_native.min20.annotated.bed | v6 | S000059 cexome6 CNV targets | deferred until n>=30 |

Format: `chr, start, end, gc, gene`  GC fraction in col4, gene symbol in col5,
sorted in ngs-bits order (chr1..chr22, chrX, chrY) via BedSort.

## Construction

All three: minimum region size 20 bp; regions with non-numeric GC dropped
(chrY PAR is N-masked in the reference, so GC is undefined there).

**legacy**  v1/v2 intervals retained whole where >=90% covered by v3/v4
(`bedtools intersect -f 0.90 -u`), rather than a plain intersection, which
would have fragmented target boundaries into ~20k sub-50bp slivers. Gene names
taken from v1/v2, gaps filled from v3/v4 (5,445 unnamed -> 129). Retains
99.25% of the v1/v2 footprint; 40 genes lost entirely, all pseudogenes,
lncRNAs, or hypervariable immune loci.

**v5 / v6**  native vendor designs, unmodified except for the size filter and
GC handling. The vendor annotation column (nested `GENE|GENE;ENST...;ClinID-...`)
was reduced to a single symbol, preferring a readable gene name over ENSG/LOC
identifiers. `backbone` labels are retained as-is; these are non-genic tiles
placed for CNV sensitivity.

## Known limitations

**Legacy has no backbone tiling.** v1-v4 were SNV designs; CNV resolution is
bounded by exon targets alone. Multi-exon and whole-gene events call reliably;
single-exon events and precise breakpoints do not. This is a property of the
assay and cannot be recovered analytically.

**Results are not comparable across groups.** Different designs, different
evidence. Do not compare cohort-level CNV frequencies between groups.

**v6 relabelled many small exon targets as backbone.** Coverage improved (larger
tiles) but labels became less specific  e.g. AGTR1's 40bp v5 target sits inside
a 600bp v6 region labelled `backbone`. Gene assignment for reporting comes from
VEP downstream, so this is cosmetic.

**v6 omits 18 regions (13.6 kb) present in v5**, almost all imprinted non-coding
RNA targets in the chr14q32 DLK1-MEG3 domain (MEG3, MEG8, SNORD114), plus
KCNQ1DN (chr11p15) and HCG22. No protein-coding gene is affected. Relevant only
if imprinting disorders are in scope.

## Rebuilding

`build_native_bed.sh` regenerates v5 and v6 from the vendor BEDs. The legacy BED
requires the v1/v3 intersect step; see the bed audit log.

Sort order matters: bedtools sweep operations need lexicographic sort, ngs-bits
needs natural chromosome order. BedSort must be the last step.
EOF

echo "wrote PROVENANCE.md"
hr
ls -la "$BEDDIR"
