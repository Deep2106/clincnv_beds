#!/usr/bin/env bash
# =============================================================================
# build_native_bed.sh - Build a ClinCNV-ready BED from a native capture design
# =============================================================================
# Produces: chr, start, end, gc, gene  (BedSort order, no n/a GC rows)
#
# Usage:
#   ./build_native_bed.sh --in  S000054_cexome5_hg38_cnvs_targets.bed \
#                         --out v5_native.min20.annotated.bed \
#                         [--min-bp 20]
#
# Env (or edit defaults below):
#   SIF_DIR   directory holding ngsbits.sif
#   FASTA     reference genome
# =============================================================================

set -euo pipefail

SIF_DIR="${SIF_DIR:-path/to/Docker_Images}"
FASTA="${FASTA:-path/to/reference/hg38.fasta}"
MIN_BP=20
IN=""; OUT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --in)     IN="$2"; shift 2 ;;
    --out)    OUT="$2"; shift 2 ;;
    --min-bp) MIN_BP="$2"; shift 2 ;;
    -h|--help) sed -n '2,18p' "$0"; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

[[ -z "$IN"  ]] && { echo "ERROR: --in required"  >&2; exit 2; }
[[ -z "$OUT" ]] && { echo "ERROR: --out required" >&2; exit 2; }
[[ -f "$IN" ]]    || { echo "ERROR: not found: $IN" >&2; exit 1; }
[[ -f "$FASTA" ]] || { echo "ERROR: not found: $FASTA" >&2; exit 1; }

NGSBITS="$SIF_DIR/ngsbits.sif"
[[ -f "$NGSBITS" ]] || { echo "ERROR: not found: $NGSBITS" >&2; exit 1; }

WORK=$(mktemp -d); trap 'rm -rf "$WORK"' EXIT
hr() { printf '%s\n' "----------------------------------------------------------"; }

echo "input   : $IN"
echo "output  : $OUT"
echo "min bp  : $MIN_BP"
hr

# --- 1. normalise: primary contigs, keep annotation, sort -------------------
grep -v '^#\|^track\|^browser' "$IN" \
  | awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|X|Y)$/ && $2 ~ /^[0-9]+$/ {
      print $1, $2, $3, ($4==""?".":$4) }' \
  | LC_ALL=C sort -k1,1 -k2,2n > "$WORK/norm.bed"
echo "1. normalised            : $(wc -l < "$WORK/norm.bed") regions"

# --- 2. reduce the annotation to a single usable symbol ---------------------
# Native v5/v6 col4 looks like:
#   backbone
#   OR4F5
#   AKT1|AKT1;ENST00000554826.2;ClinID-493153;NM_005163.2;...
#   ENSG00000236601;ENSG00000237094;LOC112268260;ENST00000450983.1;...
# Take the first token that is a plain symbol; fall back to the first token.
awk 'BEGIN{OFS="\t"}
{
  raw = $4
  if (raw == "." || raw == "" || tolower(raw) == "backbone") { print $1,$2,$3,raw; next }
  n = split(raw, tok, /[|;,]/)
  pick = ""; fallback = ""
  for (i = 1; i <= n; i++) {
    t = tok[i]
    gsub(/^[ \t]+|[ \t]+$/, "", t)
    if (t == "" || t == ".") continue
    # skip transcript / protein / clinvar / CCDS accessions outright
    if (t ~ /^(ENST|ENSP|ENSE)[0-9]/)      continue
    if (t ~ /^(NM_|NR_|XM_|XR_|NP_|XP_)/)  continue
    if (t ~ /^CCDS[0-9]/)                  continue
    if (t ~ /^ClinID-/)                    continue
    # ENSG / LOC ids are usable but a readable symbol is preferred
    if (t ~ /^ENSG[0-9]/ || t ~ /^LOC[0-9]/) { if (fallback == "") fallback = t; continue }
    pick = t; break
  }
  if (pick == "") pick = fallback
  if (pick == "") { pick = tok[1]; gsub(/^[ \t]+|[ \t]+$/, "", pick) }
  if (pick == "") pick = "."
  print $1, $2, $3, pick
}' "$WORK/norm.bed" > "$WORK/named.bed"

echo "2. annotation reduced"
echo "   backbone regions      : $(awk '$4=="backbone"' "$WORK/named.bed" | wc -l)"
echo "   unnamed (.)           : $(awk '$4=="."' "$WORK/named.bed" | wc -l)"
echo "   ENSG-only symbols     : $(awk '$4 ~ /^ENSG/' "$WORK/named.bed" | wc -l)"
echo "   distinct symbols      : $(cut -f4 "$WORK/named.bed" | grep -viE '^(backbone|\.)$' | sort -u | wc -l)"

# --- 3. merge self-overlaps, collapsing names ------------------------------
PRE=$(wc -l < "$WORK/named.bed")
bedtools merge -i "$WORK/named.bed" -c 4 -o distinct \
  | awk 'BEGIN{OFS="\t"}
    { n=$4
      if (n != ".") { sub(/^\.,/,"",n); sub(/,\.$/,"",n); gsub(/,\.,/,",",n) }
      # if a real symbol merged with backbone, prefer the symbol
      if (n ~ /,/ && n ~ /backbone/) { gsub(/^backbone,/,"",n); gsub(/,backbone$/,"",n); gsub(/,backbone,/,",",n) }
      print $1,$2,$3,n }' > "$WORK/merged.bed"
POST=$(wc -l < "$WORK/merged.bed")
echo "3. merged                : $POST regions (was $PRE; $((PRE-POST)) collapsed)"

# --- 4. minimum region size -------------------------------------------------
awk -v m="$MIN_BP" 'BEGIN{OFS="\t"} ($3-$2)>=m' "$WORK/merged.bed" > "$WORK/minbp.bed"
echo "4. >=${MIN_BP}bp filter        : $(wc -l < "$WORK/minbp.bed") regions, \
$(awk '{s+=$3-$2}END{print s}' "$WORK/minbp.bed") bp \
(dropped $(( POST - $(wc -l < "$WORK/minbp.bed") )))"

# --- 5. GC annotation -------------------------------------------------------
singularity exec -B "$WORK" -B "$(dirname "$FASTA")" "$NGSBITS" \
  BedAnnotateGC -in "$WORK/minbp.bed" -out "$WORK/gc_raw.bed" -ref "$FASTA"

GCCOL=$(awk '!/^#/ {print NF; exit}' "$WORK/gc_raw.bed")
[[ "$GCCOL" -eq 5 ]] || { echo "ERROR: expected 5 cols from BedAnnotateGC, got $GCCOL" >&2; exit 1; }
echo "5. GC annotated          : $(wc -l < "$WORK/gc_raw.bed") regions"

# --- 6. reorder to chr/start/end/gc/gene, drop non-numeric GC ---------------
awk 'BEGIN{OFS="\t"} $5 ~ /^[0-9]*\.?[0-9]+$/ {print $1,$2,$3,$5,$4}' \
  "$WORK/gc_raw.bed" > "$WORK/reordered.bed"
DROPPED=$(( $(wc -l < "$WORK/gc_raw.bed") - $(wc -l < "$WORK/reordered.bed") ))
echo "6. reordered (gc=col4)   : $(wc -l < "$WORK/reordered.bed") regions"
echo "   dropped n/a GC        : $DROPPED"
[[ "$DROPPED" -gt 0 ]] && awk '$5 !~ /^[0-9]*\.?[0-9]+$/ {print "     "$1}' "$WORK/gc_raw.bed" \
  | sort | uniq -c | sort -rn

# --- 7. BedSort (ngs-bits natural chromosome order) ------------------------
singularity exec -B "$WORK" "$NGSBITS" \
  BedSort -in "$WORK/reordered.bed" -out "$WORK/sorted.bed"

# --- 8. validate ------------------------------------------------------------
hr
BAD=$(awk 'NF!=5' "$WORK/sorted.bed" | wc -l)
[[ "$BAD" -eq 0 ]] || { echo "ERROR: $BAD rows without 5 columns" >&2; exit 1; }
BADGC=$(awk '$4 !~ /^[0-9]*\.?[0-9]+$/ || $4<0 || $4>1' "$WORK/sorted.bed" | wc -l)
[[ "$BADGC" -eq 0 ]] || { echo "ERROR: $BADGC rows with invalid GC" >&2; exit 1; }

cp "$WORK/sorted.bed" "$OUT"
echo "WROTE $OUT"
echo "  regions : $(wc -l < "$OUT")"
echo "  bases   : $(awk '{s+=$3-$2}END{print s}' "$OUT")"
echo "  md5     : $(md5sum "$OUT" | cut -d' ' -f1)"
echo
echo "chromosome order:"
cut -f1 "$OUT" | uniq | tr '\n' ' '; echo
echo
echo "first rows:"
head -3 "$OUT"
hr
echo "NEXT: verify BedCoverage accepts it (expects 6-col output):"
echo "  singularity exec -B /home/data/human_genetics -B /home/sdata \\"
echo "    $NGSBITS BedCoverage -bam <bam> -in $OUT -min_mapq 5 -decimals 4 | head -3"
