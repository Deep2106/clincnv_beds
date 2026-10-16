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
