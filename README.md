# A 600-Million-Year-Old +2 Dual-Coding Gene Family Inside Human *PRRC2A* and *PRRC2B*

Self-contained code, compact raw alignment and gnomAD v4.1 datasets, ClinVar and CBE knockout switch tables, editable Fignimble vector figures, and compiled preprint for:

> **An Uncharacterized 600-Million-Year-Old +2 Dual-Coding Gene Family Embedded Inside the Human Leaky-Scanning Regulators *PRRC2A* and *PRRC2B***  
> Pratik Koppikar (University of Cambridge)

---

## Summary of Key Findings

1. **Uncharacterized human +2 dual-coding paralog family spanning coding exons 2–6:** A frame-agnostic screen of 23,183 canonical human protein-coding transcripts (`knownCanonical`, hg38) across 100 vertebrates identifies a deeply conserved +2 overlapping reading frame inside the unstructured N-terminal "BAT2" domain (`InterPro IPR009738`) of both *PRRC2A* (chr6 MHC class III; `alt-PRRC2A`, 130 aa, CDS nt 254–644) and *PRRC2B* (chr9; `alt-PRRC2B`, 142 aa, CDS nt 257–683; raw rank 2 genome-wide). Both proteins share 59% amino-acid identity across their aligned 129-residue core (76/129; 53% overall) across three conserved blocks (`SAAG[I/V]AA[T/A]AGF[T/A]`, `KG[I/L]KPTV[T/I][I/L]LSRGIS`, and `WR[P/A]GQG...VVWARTKPPP`). Although HAVANA/GENCODE quietly added unnamed stubs on May 25, 2022 (`ENSG00000289282` / UniProt `A0A8I5QKQ9`, `PE=1` via PeptideAtlas; `ENSG00000288841` / `A0A8I5KXM2`, `PE=4`), both lack HGNC symbols, NCBI Gene/RefSeq/GTEx entries, Pfam/InterPro domains, and any published characterization, with Ensembl Compara listing only human and mouse.
2. **600-My deuterostome origin and paralog-specific loss in *PRRC2C*:** Absent in protostomes (*Drosophila* `nocte`, 4 stops in +2), the +2 frame opens at the base of deuterostomes (~580 Mya; 0 stops in acorn worm `XP_077866102.1` and sea urchin `XP_030852947.1`), locks its alanine-rich and `VVWARTKPPP` motifs through amphioxus (`XP_035683835.1`) and sea lamprey (`XP_032802661.2`), and survives early vertebrate 2R whole-genome duplications in *PRRC2A* and *PRRC2B* while pseudogenizing in *PRRC2C* prior to jawed vertebrates (>465 Mya; 13 stops in elephant shark `prrc2c` `XM_042335521.1` and 10 stops in human `PRRC2C`).
3. **Orthogonal vertebrate purifying selection against built-in paralog negative control *PRRC2C*:** Across Panel-32 vertebrates, strict non-stop-risk fourfold degenerate (4D) synonymous identity surges inside the +2 overlap (*PRRC2A* `0.874` vs. `0.613`, $\Delta = +0.262$, 22/22 species, exact one-sided paired Wilcoxon $P = 2.38 \times 10^{-7}$; *PRRC2B* `0.868` vs. `0.570`, $\Delta = +0.298$, 26/26 species, $P = 1.49 \times 10^{-8}$; replicating in 8/8 major vertebrate clades, both $P = 3.91 \times 10^{-3}$; both $P \le 3.74 \times 10^{-6}$ vs. *PRRC2C* $\Delta = +0.123$). Only Frame +2 collapses to near-zero stops (`0.05` and `0.08` vs. `6.15` and `5.97` expected under a species-specific within-gene flanking-codon synonymous null, $P \le 2.38 \times 10^{-7}$; within-gene max $P \le 1.51 \times 10^{-5}$ vs. all 4 other frames [1.00–7.52 stops]; corroborated by a within-window codon-bag shuffle on hg38, $P \le 4.00 \times 10^{-3}$; `94/98` *PRRC2B* and `74/84` *PRRC2A* 100-way vertebrates stop-free vs. `0/85` in *PRRC2C*). Where Frame-0 amino acids diverge from human, the overlapping +2 residue(s) directly touched by the mutated nucleotide(s) stay identical at `28.8%` (*PRRC2A*) and `53.5%` (*PRRC2B*) vs. `8.3%` (*PRRC2C*; $P = 3.14 \times 10^{-5}$ and $3.39 \times 10^{-7}$).
4. **Active purifying selection across 807,162 human exomes and genomes (gnomAD v4.1, QC-PASS):** Holding Frame 0 strictly synonymous (`F0 = syn`), **0/35** non-CpG Frame-0-silent / Frame-+2-nonsense mutations are observed across pooled *PRRC2A*+*PRRC2B* (`0.0%`, `AC = 0`) vs. **94/533** non-stop sites (`17.6%`, `AC = 62,444`; exact Fisher one-sided $P = 1.43 \times 10^{-3}$, two-sided $P = 1.96 \times 10^{-3}$), alongside hypermutable CpG depletion (`2/7`, `AC = 3` vs. `25/28`, `AC = 476`; $P = 3.04 \times 10^{-3}$; three-stratum exact hypergeometric CMH common $\text{OR} = 0.019$, exact conditional MLE $\text{OR} = 0.041$ [95% CI `0.003–0.261`], exact one-sided $P = 9.06 \times 10^{-6}$, exact two-sided $P = 1.28 \times 10^{-5}$), replicating independently in *PRRC2A* ($P = 1.15 \times 10^{-3}$) and *PRRC2B* ($P = 1.10 \times 10^{-2}$) while disappearing in *PRRC2C* (`3/17 = 17.6%`, $\text{OR} = 0.763$, exact $P = 0.490$). Among Frame-0-missense sites (`F0 = mis`), observation rates follow a monotonic purifying gradient (`+2 SYN > +2 MIS > +2 STOP`; exact CMH `+2 MIS` $\text{OR} = 0.645$, $P = 1.47 \times 10^{-3}$; `+2 STOP` $\text{OR} = 0.201$, $P = 5.34 \times 10^{-5}$; null in *PRRC2C*).
5. **Constitutive splicing, Kozak asymmetry, ancestral +2 uoORF reinitiation architecture, and genomically audited knockout switches:** GTEx v2 Snaptron junctions (~19,000 samples) across all 5 symmetric canonical introns spanning coding exons 1–6 rule out alternative splicing (`0.52%` and `1.05%` top internal alternative junction vs. mean canonical + alt reads). Across 100 vertebrates, Frame-0 preserves a 100%-conserved suboptimal `+4 T` Kozak base (`98/98`, `99/99`, `100/100`; adequate `G/ATG/T` in all three human paralogs) while +2 preserves 100%-optimal `[A/G]NNATGG` (`98/98` in *PRRC2A*, `97/97` ATG-bearing [97/98 covered] in *PRRC2B*, relaxed to `60/70` [85.7%] in *PRRC2C*). GENCODE Phase I Ribo-seq identifies a replicated ancestral +2 uoORF (`c6riboseqorf41`, `c9riboseqorf98`, `c1riboseqorf224`) at `-16 nt` terminating exactly **183 nt** upstream of the +2 start codon in both *PRRC2A* and *PRRC2B*—matching the canonical *ATF4* reinitiation window. We re-annotate 23 ClinVar variants altering Frame +2 and provide genomically audited single-nucleotide knockout switches—including bystander-free `SpCas9-NG`/`SpRY`-CBE `Q4*` switches (`PRRC2A c.264C>T`, `PRRC2B c.267C>T`), bystander-free canonical `SpCas9-NGG`-CBE `Q96*`/`Q97*` switches, and prime-editing targets where adjacent `-1/-2` cytosines carry Frame-0 missense risk.

---

## Repository Layout

```text
prrc2-plus2-dualcoding/
├── preprint/
│   ├── prrc2_plus2_preprint.tex           # Preprint LaTeX manuscript (pure black Arial typography)
│   ├── prrc2_plus2_preprint.pdf           # Compiled 15-page vector PDF
│   └── figures/                           # Native Fignimble (pachterlab/Fignimble) project
│       ├── style.py                       # Shared Nature styling and matte semantic palette
│       ├── graphical_abstract.py          # Graphical Abstract
│       ├── fig1.py ... fig4.py            # Main Figures 1-4 (183 mm double-column layout)
│       ├── panels/                        # 16 modular, data-driven Matplotlib panel scripts
│       ├── rendered/                      # Exported vector PDFs
│       └── build_all.sh                   # One-command preflight check + export script
├── data/
│   ├── prrc2_100way_exons.fa.gz           # Exact 100-vertebrate UCSC exon alignments (PRRC2A/B/C)
│   ├── gnomad_v4_prrc2_variants.json.gz   # Exact gnomAD v4.1 region SNVs (with QC filters) & pLI/LOEUF
│   ├── gtex_v2_snaptron_junctions.json    # Exact GTEx v2 Snaptron splice junctions (coding exons 1-6)
│   ├── metazoan_prrc2_cds.fa              # Exact 7 non-human NCBI RefSeq metazoan PRRC2/nocte CDSs
│   ├── hg38_prrc2_switch_flanks.json      # Exact hg38 genomic flanks for Kozak & CBE PAM/bystander audit
│   ├── clinvar_prrc2_dual_summary.json    # ClinVar variants across PRRC2A/B coding exons 2-6
│   ├── gencode_riboseq_prrc2_uoorfs.bed   # GENCODE Phase I Ribo-seq uoORF annotations (hg38)
│   ├── track1_screen_summary.json         # Genome-wide 23,183-transcript 100-way screen summary
│   ├── uniprot_and_compara/               # UniProt, Ensembl Compara, MMseqs2, and fish audit JSONs
│   └── master_preprint_data.json          # Unified verified dataset regenerated by script 01
├── scripts/
│   ├── 01_reproduce_all_analyses.py       # 100% self-contained pipeline (recomputes all stats in ~5 s)
│   ├── 02_genome_wide_100way_scan.py      # Upstream 23,183-transcript UCSC 100-way scanner
│   └── 03_trace_metazoan_origin.py        # Live NCBI RefSeq verification script for all 7 metazoan lineages
├── CITATION.cff
└── LICENSE
```

---

## Reproducing All Analyses, Figures, and the Preprint

### 1. Re-run the full self-contained analysis pipeline (zero external files required)

```bash
python3 scripts/01_reproduce_all_analyses.py
```

*(Requires Python 3.10+ with `numpy` and `scipy`; regenerates `data/master_preprint_data.json` directly from the compact raw subsets in `data/` in ~5 seconds.)*

### 2. Re-render or visually inspect figures with Fignimble

```bash
uv tool install "git+https://github.com/pachterlab/Fignimble.git"
bash preprint/figures/build_all.sh
```

### 3. Compile the preprint PDF

```bash
cd preprint
tectonic prrc2_plus2_preprint.tex
```

---

## License and Data Provenance

- Code, figures, and manuscript text: **MIT License** (`LICENSE`).
- Upstream public reference datasets used to generate the compact subsets in `data/`: UCSC 100-way `knownCanonical.exonNuc.fa.gz` (hg38), gnomAD v4.1 GraphQL API (`CC0`), GENCODE Phase I Ribo-seq ORF catalog (`GENCODE_Riboseq_ORFs.bed`), NCBI ClinVar and RefSeq, UniProtKB (`CC-BY-4.0`), Ensembl Compara, and GTEx v2 Snaptron.
