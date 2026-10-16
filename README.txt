ClinCNV BEDs - cohort split by capture kit
==========================================
Location : path/to/reference/EXOME_XtraXXX/clincnv_beds/
Built    : 2026-09-11 (build_native_bed.sh, PROVENANCE.md, checksums.txt in same folder)


1. WHY THE COHORT IS SPLIT INTO 3 GROUPS
----------------------------------------
ClinCNV calls CNVs by comparing read depth across a reference cohort. Samples
captured with different kits cannot be pooled: a region missing from one kit's
design has near-zero depth and is called as a homozygous deletion, not as
"not tested".

The six kit versions are really three distinct designs:
  - v1 = v2 (identical md5), v3 = v4 (identical md5): SNV targets only, no CNV BED, no backbone
  - v5 (S000054, cexome5): dedicated CNV BED, exons + backbone tiling
  - v6 (S000059, cexome6): dedicated CNV BED, exons + more backbone tiling (+15.9 Mb vs v5)

v5 and v6 were not merged. A v5 intersect v6 consensus keeps 99.89% of v5 but
removes about 19% of v6's design (mostly backbone tiling).

As per date 14/10/2026:

  Group   Kits                                         Samples  Status
  legacy  v1-v4 (S000032, S000042, S000048, S000052)   71       ready to run
  v5      S000054                                      39       ready to run
  v6      S000059                                      20       DEFERRED until ~30 samples
                                                                (too few for a stable reference cohort)

Within each group, trios and singletons/duos run as separate ClinCNV pipelines.
The coverage merge (--normal) uses ALL samples in the group, whatever their
family structure.


2. BED FILES AND HOW THEY WERE BUILT
------------------------------------
  legacy_consensus.min20.annotated.bed   207,716 regions   ~36.2 Mb
    - Source: v1/v2 and v3/v4 SNV target BEDs
    - bedtools intersect -f 0.90 -u: keeps a v1/v2 interval whole if >=90%
      of it overlaps v3/v4, so targets are not fragmented
    - Gene names missing (".") in v1/v2 filled from v3/v4

  v5_native.min20.annotated.bed          265,192 regions   66,987,205 bp
  v6_native.min20.annotated.bed          272,258 regions   82,862,969 bp
    - Source: each kit's own vendor CNV BED (native, no intersect)
    - Vendor CNV BEDs are non-overlapping (0 regions collapsed on merge)

Common processing (all three):
  a. Drop regions < 20 bp (params.min_target_bp = 20)
     - <20 bp regions average ~10-12 bp: one or two read-depth observations,
       dominated by Poisson noise and edge effects, with unreliable GC.
       They inflate ClinCNV's per-region variance.
     - Cost: 0.06-0.17% of footprint (v1: 6,211, v5: 6,364, v6: 4,018 regions)
     - 50 bp rejected: it costs 0.6-1.7% of footprint, and 20-49 bp
       regions are still usable.
  b. BedAnnotateGC + BedAnnotateGenes (ngs-bits), using
     hg38_v0_Homo_sapiens_assembly38.fasta
  c. Reorder to chr / start / end / gc / gene (ClinCNV expects GC in col4;
     the vendor col4 annotation is discarded)
  d. Drop regions with n/a GC (N-masked chrY PAR)
  e. BedSort chromosome order (chr1..chr22, chrX, chrY)


3. ClinCNV SETTINGS (all groups)
--------------------------------
  --scoreG 50 --lengthG 0 --minimumNumOfElemsInCluster 10
  --maxNumGermCNVs 2000 --maxNumIter 5
  --par "chrX:10001-2781479;chrX:155701383-156030895"
    chrY PAR is N-masked, so all PAR reads map to chrX. Without --par,
    every male sample gets spurious duplications in these regions.


4. CAVEATS
----------
  - Calls are NOT directly comparable across groups (different callable regions).
  - A gene outside a group's BED is "no call", NOT "no CNV". Reports must say so.
  - Changing a sample's BED changes its calls. Publish the BED used with every
    run, and version output paths (e.g. clincnv/v5_v1/).
  - Some <20 bp regions are deliberate small-exon or hotspot targets that are
    no longer callable.
  - OPEN: check backbone depth on a real v5 (S000054) BAM before production.
    If it is ~0, the backbone was not physically captured and v5/v6 must fall
    back to exon-only targets.
  - v6: re-evaluate when it reaches ~30 samples, then run as its own cohort.
    This is a deliberate re-analysis, not a silent change.
