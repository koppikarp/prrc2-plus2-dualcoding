#!/usr/bin/env python3
"""Deterministically build Supplementary Information LaTeX/PDF for the
PRRC2A/B +2 dual-coding preprint directly from verified repository datasets.
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
PREPRINT = ROOT / "preprint"


def esc(s):
    if s is None:
        return "--"
    s = str(s)
    for k, v in [("\\", r"\textbackslash{}"), ("_", r"\_"), ("%", r"\%"),
                 ("&", r"\&"), ("#", r"\#"), ("$", r"\$"), (">", r"$>$"), ("<", r"$<$")]:
        s = s.replace(k, v)
    return s


def fmt_p(p):
    if p is None:
        return "--"
    if p < 1e-3:
        m, e = f"{p:.2e}".split("e")
        return f"${m} \\times 10^{{{int(e)}}}$"
    return f"{p:.4f}"


def build_supplement_tex():
    screen = json.loads((DATA / "track1_screen_summary.json").read_text())
    pm = json.loads((DATA / "master_preprint_data.json").read_text())

    # Table S1: Top 35 from 23,183-gene screen
    t1_rows = []
    for rank, h in enumerate(screen["top35"], 1):
        gene = h["gene"]
        if h["tx_id"] == "ENST00000357304.8":
            gene = "PRRC2B (ENST00000357304.8)"
        elif h["tx_id"] == "ENST00000454306.6":
            gene = "PRRC2A (ENST00000454306.6)"
        is_prrc = "PRRC2" in gene
        gfmt = f"\\textbf{{{esc(gene)}}}" if is_prrc else esc(gene)
        flag = "Single-exon" if h["single_exon"] else f"{h['n_exons_spanned']} exons"
        t1_rows.append(
            f"{rank} & {gfmt} & {esc(h['frame'])} & {h['alt_len_aa']} & "
            f"{h['id4d_in']:.3f} & {h['id4d_out']:.3f} & $+{h['d4d']:.3f}$ & "
            f"{h['sp_stopfree']}/{h['sp_covered']} & {esc(flag)} & {h['score']:.1f} \\\\"
        )
    t1_lines = "\n".join(t1_rows)

    # Table S2a: 10 metazoan lineages spanning 680 Myr
    t2a_lines = "\n".join(
        f"{esc(r['display_label'])} & {esc(r['clade'])} & {r['mya']} & "
        f"\\texttt{{{esc(r['accession'])}}} & {r['stops_p2_core']} & {r['stops_p1_core']} & "
        f"{'Intact' if r['plus2_intact'] else 'Lost/Absent'} & \\texttt{{{esc(r['motif_11aa'])}}} \\\\"
        for r in pm["evolutionary_origin"]
    )

    # Table S2b: Comparative vertebrate selection audit
    vs = pm["vertebrate_selection"]
    t2b_rows = []
    for g in ["PRRC2A", "PRRC2B", "PRRC2C"]:
        v = vs[g]
        cl = v["nsr4d_clade_collapsed"]
        f2 = v["five_frame_stops_panel32"]["+2"]
        cb = v["hg38_codon_bag_null"]["+2"]
        t2b_rows.append(
            f"\\gene{{{g}}} & {v['n_panel32_4d']} & {v['nsr4d_in_mean']:.3f} / {v['nsr4d_out_mean']:.3f} "
            f"($+{v['nsr4d_delta']:.3f}$) & {fmt_p(v['nsr4d_p_one'])} & "
            f"$+{cl['delta']:.3f}$ ({fmt_p(cl['p_one'])}) & "
            f"{f2['obs_mean']:.2f} vs.\\ {f2['exp_mean']:.2f} ({fmt_p(f2['p_less_vs_null'])}) & "
            f"{cb['obs']} vs.\\ {cb['exp_mean']:.2f} ({fmt_p(cb['emp_p_le'])}) & "
            f"{v['all100_stopfree_n']}/{v['all100_covered_n']} & "
            f"{v['p2_aa_pres_at_f0_div_mean']*100:.1f}\\% \\\\"
        )
    t2b_lines = "\n".join(t2b_rows)

    # Table S3a: gnomAD v4.1 per-stratum counts
    gn = pm["gnomad_v4"]
    cohorts = [
        ("PRRC2A", "PRRC2A"),
        ("PRRC2B", "PRRC2B"),
        ("POOLED_AB", "Pooled A+B"),
        ("POOLED_AB_NO_EDGE2BP", "Pooled A+B (no $\\pm 2$ bp)"),
        ("PRRC2C_NEGCTRL", "PRRC2C (neg.\\ ctrl)"),
    ]
    t3a_rows = []
    for key, lbl in cohorts:
        st = gn[key]["strata"]
        for strat, slbl in [("CpG_Ti", "CpG Ti"), ("nonCpG_Ti", "non-CpG Ti"), ("Tv", "Tv")]:
            s = st[strat]
            syn_ns_obs = s["F0_syn__p2_syn"]["obs"] + s["F0_syn__p2_mis"]["obs"]
            syn_ns_tot = s["F0_syn__p2_syn"]["tot"] + s["F0_syn__p2_mis"]["tot"]
            syn_st_obs = s["F0_syn__p2_stop_gained"]["obs"]
            syn_st_tot = s["F0_syn__p2_stop_gained"]["tot"]
            syn_st_ac = s["F0_syn__p2_stop_gained"]["sumAC"]
            mis_syn = s["F0_mis__p2_syn"]
            mis_mis = s["F0_mis__p2_mis"]
            mis_st = s["F0_mis__p2_stop_gained"]
            t3a_rows.append(
                f"{lbl} & {slbl} & "
                f"{syn_ns_obs}/{syn_ns_tot} & {syn_st_obs}/{syn_st_tot} (AC={syn_st_ac}) & "
                f"{mis_syn['obs']}/{mis_syn['tot']} & {mis_mis['obs']}/{mis_mis['tot']} & "
                f"{mis_st['obs']}/{mis_st['tot']} \\\\"
            )
    t3a_lines = "\n".join(t3a_rows)

    # Table S3b: Exact CMH & exact conditional MLE OR summary
    t3b_rows = []
    for key, lbl in cohorts:
        for test_key, tlbl in [
            ("cmh_F0syn_p2STOP_vs_nonSTOP", "F0-syn: +2 STOP vs.\\ non-STOP"),
            ("cmh_F0mis_p2MIS_vs_SYN", "F0-mis: +2 MIS vs.\\ +2 SYN"),
            ("cmh_F0mis_p2STOP_vs_SYN", "F0-mis: +2 STOP vs.\\ +2 SYN"),
        ]:
            c = gn[key][test_key]
            t3b_rows.append(
                f"{lbl} & {tlbl} & {c['common_or']:.3f} & {c['exact_mle_or']:.3f} & "
                f"[{c['exact_ci_lo']:.3f}, {c['exact_ci_hi']:.3f}] & "
                f"{fmt_p(c['p_exact_one_sided_less'])} & {fmt_p(c['p_exact_two_sided'])} \\\\"
            )
    t3b_lines = "\n".join(t3b_rows)

    # Table S3c: GTEx v2 Snaptron splice junctions
    gtex = pm["translational_architecture"]["gtex_v2_snaptron"]
    t3c_rows = []
    for g in ["PRRC2A", "PRRC2B"]:
        gd = gtex[g]
        for intr in gd["canonical_introns"]:
            t3c_rows.append(
                f"\\gene{{{g}}} & Canonical Intron {esc(intr['intron'])} & "
                f"\\texttt{{{esc(intr['chrom'])}:{intr['start']:,}-{intr['end']:,}}} & "
                f"{intr['reads']:,} & {intr['samples']:,} \\\\"
            )
        alt_in = gd["top_alt_internal_junction"]
        alt_ex1 = gd["top_alt_skip_exon1_junction"]
        t3c_rows.append(
            f"\\gene{{{g}}} & Top internal alt.\\ junction & "
            f"\\texttt{{{esc(alt_in['chrom'])}:{alt_in['start']:,}-{alt_in['end']:,}}} & "
            f"{alt_in['reads']:,} & {alt_in['samples']:,} \\\\"
        )
        t3c_rows.append(
            f"\\gene{{{g}}} & Top Exon-1-skipping junction & "
            f"\\texttt{{{esc(alt_ex1['chrom'])}:{alt_ex1['start']:,}-{alt_ex1['end']:,}}} & "
            f"{alt_ex1['reads']:,} & {alt_ex1['samples']:,} \\\\"
        )
    t3c_lines = "\n".join(t3c_rows)

    # Table S4: All 45 ClinVar variants
    cv = sorted(pm["clinvar"], key=lambda x: (x["gene"], x["p2_aa_pos"], x["variant_id"]))
    assert len(cv) == 45, f"Expected 45 ClinVar SNVs, got {len(cv)}"
    t4_rows = []
    for r in cv:
        alters = r["p2_eff"] != "SYN"
        p2_str = f"\\textbf{{\\texttt{{{esc(r['p2_change'])}}} ({esc(r['p2_eff'])})}}" if alters else f"\\texttt{{{esc(r['p2_change'])}}} ({esc(r['p2_eff'])})"
        t4_rows.append(
            f"\\gene{{{r['gene']}}} & \\texttt{{{esc(r['variant_id'])}}} & \\texttt{{{esc(r['hgvs_c'])}}} & "
            f"{esc(r['clin_sig'])} & \\texttt{{{esc(r['f0_change'])}}} ({esc(r['f0_eff'])}) & "
            f"{p2_str} & {r['p2_aa_pos']} \\\\"
        )
    t4_lines = "\n".join(t4_rows)

    # Table S5: All 42 single-nucleotide knockout switches
    sw_all = pm["experimental_switches"]["PRRC2A"] + pm["experimental_switches"]["PRRC2B"]
    assert len(sw_all) == 42, f"Expected 42 switches, got {len(sw_all)}"
    t5_rows = []
    for s in sw_all:
        pams = []
        for p in s.get("ngg_pams_pos4_8", []):
            pams.append(f"NGG:{p['pam']}(p{p['pos']})")
        for p in s.get("ng_pams_pos4_8", []):
            if not any(x["pam"] == p["pam"] and x["pos"] == p["pos"] for x in s.get("ngg_pams_pos4_8", [])):
                pams.append(f"NG:{p['pam']}(p{p['pos']})")
        pam_str = ", ".join(pams[:2]) if pams else "--"
        bys = s.get("bystanders_minus3_plus3", [])
        if not s["is_ct_transition"]:
            bys_str = "N/A"
        elif not bys:
            bys_str = "\\textbf{None}"
        else:
            bys_str = ", ".join(f"{b['offset']:+d}C({b['f0_eff']})" for b in bys[:3])
        ed_short = s["editor"].replace("Prime Editing / HDR", "PE / HDR").replace("Prime Editing (PE) / narrow-window CBE", "PE / narrow CBE").replace("SpCas9-NG / SpRY-CBE or PE", "NG/SpRY-CBE or PE")
        t5_rows.append(
            f"\\gene{{{s['gene']}}} & \\texttt{{c.{s['cds_nt']}{esc(s['sub'])}}} & "
            f"\\texttt{{{esc(s['chrom'])}:{s['hg38_pos']:,}}} & "
            f"Ex{s['exon']} & \\texttt{{{esc(s['f0_aa'])}}} ({s['f0_usage_ref']:.1f}$\\rightarrow${s['f0_usage_alt']:.1f}) & "
            f"\\texttt{{{esc(s['p2_aa'])}}} & {esc(pam_str)} & {bys_str} & {esc(ed_short)} \\\\"
        )
    t5_lines = "\n".join(t5_rows)

    tex = rf"""\documentclass[10pt,letterpaper]{{article}}

\usepackage[letterpaper,left=0.60in,right=0.60in,top=0.72in,bottom=0.72in,headheight=19pt]{{geometry}}
\usepackage{{fontspec}}
\usepackage{{amsmath,amssymb}}
\usepackage{{booktabs}}
\usepackage{{longtable}}
\usepackage{{array}}
\usepackage{{xcolor}}
\usepackage{{titlesec}}
\usepackage{{fancyhdr}}
\usepackage{{lastpage}}
\usepackage[hidelinks]{{hyperref}}
\urlstyle{{same}}

\setmainfont{{Arial}}[
  Scale=0.92,
  BoldFont={{Arial Bold}},
  ItalicFont={{Arial Italic}},
  BoldItalicFont={{Arial Bold Italic}}
]
\setsansfont{{Arial}}[Scale=0.92]
\setmonofont{{Arial}}[Scale=0.92]

\linespread{{1.08}}
\setlength{{\parindent}}{{0pt}}
\setlength{{\parskip}}{{4pt}}
\definecolor{{RuleGray}}{{HTML}}{{BDBDBD}}
\color{{black}}

\titleformat{{\section}}
  {{\fontsize{{11.5}}{{14}}\selectfont\bfseries\color{{black}}}}
  {{}}{{0pt}}{{}}
  [\vspace{{-3pt}}\textcolor{{RuleGray}}{{\rule{{\textwidth}}{{0.5pt}}}}]
\titlespacing*{{\section}}{{0pt}}{{12pt}}{{6pt}}

\pagestyle{{fancy}}
\fancyhf{{}}
\renewcommand{{\headrulewidth}}{{0.4pt}}
\renewcommand{{\headrule}}{{\hbox to\headwidth{{\color{{RuleGray}}\leaders\hrule height \headrulewidth\hfill}}}}
\fancyhead[L]{{\footnotesize\color{{black}}Koppikar (2026) \textbar\ Supplementary Information: 600-My +2 dual-coding gene family inside human \textit{{PRRC2A}} and \textit{{PRRC2B}}}}
\fancyhead[R]{{\footnotesize\color{{black}}Supplementary Information}}
\fancyfoot[C]{{\footnotesize\color{{black}}S\thepage\ / S\pageref{{LastPage}}}}

\newcommand{{\gene}}[1]{{\textit{{#1}}}}
\newcommand{{\spname}}[1]{{\textit{{#1}}}}

\begin{{document}}

\noindent{{\fontsize{{13.5}}{{17}}\selectfont\bfseries\color{{black}}\raggedright
Supplementary Information\\[3pt]
An Uncharacterized 600-Million-Year-Old +2 Dual-Coding Gene Family Embedded Inside the Human Leaky-Scanning Regulators \gene{{PRRC2A}} and \gene{{PRRC2B}}\par}}
\vspace{{4pt}}
\noindent{{\small Pratik Koppikar\textsuperscript{{1,*}} \;\textbar\; \textsuperscript{{1}}University of Cambridge, Cambridge, United Kingdom \;\textbar\; \href{{mailto:psk34@cam.ac.uk}}{{psk34@cam.ac.uk}}\\
Repository and raw datasets: \url{{https://github.com/koppikarp/prrc2-plus2-dualcoding}}\par}}
\vspace{{3pt}}
\noindent{{\color{{RuleGray}}\rule{{\textwidth}}{{0.7pt}}}}

\section*{{Contents}}
\begin{{itemize}}
  \item \textbf{{Supplementary Table S1:}} Top 35 ranked candidates from the genome-wide 23,183-gene 100-vertebrate dual-coding screen.
  \item \textbf{{Supplementary Table S2:}} Deep metazoan evolutionary origin (10 animal lineages, 680 Myr) and complete Panel-32 / 100-way vertebrate comparative selection statistics.
  \item \textbf{{Supplementary Table S3:}} Complete gnomAD v4.1 (807,162 exomes/genomes, QC-PASS) mutation-class-stratified contingency tables, exact CMH/MLE odds ratios, and GTEx v2 Snaptron splice-junction read counts ($\sim 19,000$ samples).
  \item \textbf{{Supplementary Table S4:}} Complete dual-frame re-annotation of all 45 ClinVar coding SNVs in \gene{{PRRC2A}} and \gene{{PRRC2B}} coding exons 2--6.
  \item \textbf{{Supplementary Table S5:}} Complete catalog of all 42 genomically audited single-nucleotide silent-Frame-0 / nonsense-Frame-+2 knockout switches across \gene{{PRRC2A}} ($n=22$) and \gene{{PRRC2B}} ($n=20$).
\end{{itemize}}

\newpage
\section*{{Supplementary Table S1: Top 35 Hits From the 23,183-Gene 100-Vertebrate Screen}}
\noindent{{\footnotesize Raw uncurated top 35 hits across 23,183 canonical hg38 transcripts (3,904 deduplicated raw hits before filtering out single-exon ORFs, readthrough transcripts, and previously characterized overlapping genes). Multi-exon \gene{{PRRC2B}} (rank 2) and \gene{{PRRC2A}} (rank 29) are highlighted in bold.\par}}
\vspace{{3pt}}
{{\scriptsize\setlength{{\tabcolsep}}{{3.2pt}}\renewcommand{{\arraystretch}}{{1.08}}
\begin{{longtable}}{{@{{}}rlcrrrrrlr@{{}}}}
\toprule
\textbf{{Rank}} & \textbf{{Gene / Transcript}} & \textbf{{Frame}} & \textbf{{Alt aa}} & \textbf{{4D in}} & \textbf{{4D out}} & \textbf{{$\Delta$4D}} & \textbf{{Stop-free}} & \textbf{{Structure}} & \textbf{{Score}} \\
\midrule
\endfirsthead
\toprule
\textbf{{Rank}} & \textbf{{Gene / Transcript}} & \textbf{{Frame}} & \textbf{{Alt aa}} & \textbf{{4D in}} & \textbf{{4D out}} & \textbf{{$\Delta$4D}} & \textbf{{Stop-free}} & \textbf{{Structure}} & \textbf{{Score}} \\
\midrule
\endhead
\bottomrule
\endfoot
{t1_lines}
\end{{longtable}}
}}

\newpage
\section*{{Supplementary Table S2: Metazoan Origin and Vertebrate Comparative Selection Audit}}
\noindent{{\footnotesize\textbf{{Table S2a.}} Deep evolutionary origin and paralog-specific loss of \gene{{alt-PRRC2}} across 10 representative animal lineages spanning 680 million years.\par}}
\vspace{{3pt}}
{{\scriptsize\setlength{{\tabcolsep}}{{3.0pt}}\renewcommand{{\arraystretch}}{{1.10}}
\begin{{longtable}}{{@{{}}llrlllcl@{{}}}}
\toprule
\textbf{{Taxon / Paralog}} & \textbf{{Clade}} & \textbf{{Mya}} & \textbf{{Accession}} & \textbf{{+2 Stops}} & \textbf{{+1 Stops}} & \textbf{{+2 Status}} & \textbf{{Core 11-aa Motif}} \\
\midrule
\endhead
\bottomrule
\endfoot
{t2a_lines}
\end{{longtable}}
}}

\vspace{{8pt}}
\noindent{{\footnotesize\textbf{{Table S2b.}} Complete Panel-32 vertebrate comparative selection statistics comparing \gene{{PRRC2A}} and \gene{{PRRC2B}} against the pseudogenized paralog control \gene{{PRRC2C}} (strict non-stop-risk 4D identity, 8-clade collapse, 5-frame stop purging vs.\ flanking-codon null and hg38 within-window codon-bag null, 100-way stop-free census, and +2 amino-acid preservation at Frame-0-divergent sites).\par}}
\vspace{{3pt}}
{{\scriptsize\setlength{{\tabcolsep}}{{2.6pt}}\renewcommand{{\arraystretch}}{{1.12}}
\begin{{longtable}}{{@{{}}lclllllll@{{}}}}
\toprule
\textbf{{Gene}} & \textbf{{$n_{32}$}} & \textbf{{NSR-4D in / out ($\Delta$)}} & \textbf{{1-s $P$}} & \textbf{{8-Clade $\Delta$ ($P$)}} & \textbf{{+2 Stops vs.\ null ($P$)}} & \textbf{{Codon-bag ($P$)}} & \textbf{{100-way}} & \textbf{{+2 AA pres.}} \\
\midrule
\endhead
\bottomrule
\endfoot
{t2b_lines}
\end{{longtable}}
}}

\newpage
\section*{{Supplementary Table S3: gnomAD v4.1 Constraint and GTEx v2 Splice-Junction Audit}}
\noindent{{\footnotesize\textbf{{Table S3a.}} Mutation-class-stratified observed/possible QC-PASS SNV counts in gnomAD v4.1 (807,162 exomes and genomes) holding Frame-0 consequence strictly constant.\par}}
\vspace{{3pt}}
{{\scriptsize\setlength{{\tabcolsep}}{{3.5pt}}\renewcommand{{\arraystretch}}{{1.08}}
\begin{{longtable}}{{@{{}}lllllll@{{}}}}
\toprule
\textbf{{Cohort}} & \textbf{{Stratum}} & \textbf{{F0-syn / +2 non-STOP}} & \textbf{{F0-syn / +2 STOP}} & \textbf{{F0-mis / +2 SYN}} & \textbf{{F0-mis / +2 MIS}} & \textbf{{F0-mis / +2 STOP}} \\
\midrule
\endfirsthead
\toprule
\textbf{{Cohort}} & \textbf{{Stratum}} & \textbf{{F0-syn / +2 non-STOP}} & \textbf{{F0-syn / +2 STOP}} & \textbf{{F0-mis / +2 SYN}} & \textbf{{F0-mis / +2 MIS}} & \textbf{{F0-mis / +2 STOP}} \\
\midrule
\endhead
\bottomrule
\endfoot
{t3a_lines}
\end{{longtable}}
}}

\vspace{{6pt}}
\noindent{{\footnotesize\textbf{{Table S3b.}} Three-stratum exact hypergeometric Cochran--Mantel--Haenszel (CMH) tests and exact conditional MLE odds ratios with exact conditional 95\% confidence intervals (zero continuity hacks).\par}}
\vspace{{3pt}}
{{\scriptsize\setlength{{\tabcolsep}}{{3.2pt}}\renewcommand{{\arraystretch}}{{1.08}}
\begin{{longtable}}{{@{{}}llrrrrr@{{}}}}
\toprule
\textbf{{Cohort}} & \textbf{{Contrast}} & \textbf{{MH OR}} & \textbf{{Exact MLE OR}} & \textbf{{Exact 95\% CI}} & \textbf{{Exact 1-s $P$}} & \textbf{{Exact 2-s $P$}} \\
\midrule
\endfirsthead
\toprule
\textbf{{Cohort}} & \textbf{{Contrast}} & \textbf{{MH OR}} & \textbf{{Exact MLE OR}} & \textbf{{Exact 95\% CI}} & \textbf{{Exact 1-s $P$}} & \textbf{{Exact 2-s $P$}} \\
\midrule
\endhead
\bottomrule
\endfoot
{t3b_lines}
\end{{longtable}}
}}

\vspace{{6pt}}
\noindent{{\footnotesize\textbf{{Table S3c.}} GTEx v2 Snaptron splice-junction read support across $\sim 19,000$ human RNA-seq samples for all five canonical coding Exon 1--6 introns and highest-abundance alternative junctions.\par}}
\vspace{{3pt}}
{{\scriptsize\setlength{{\tabcolsep}}{{4.0pt}}\renewcommand{{\arraystretch}}{{1.08}}
\begin{{longtable}}{{@{{}}lllrr@{{}}}}
\toprule
\textbf{{Gene}} & \textbf{{Junction Type}} & \textbf{{hg38 Intron Coordinates}} & \textbf{{Total Reads}} & \textbf{{Samples}} \\
\midrule
\endhead
\bottomrule
\endfoot
{t3c_lines}
\end{{longtable}}
}}

\newpage
\section*{{Supplementary Table S4: Complete Re-Annotation of All 45 ClinVar SNVs in Exons 2--6}}
\noindent{{\footnotesize All 45 ClinVar coding SNVs falling within the dual-coding spans of human \gene{{PRRC2A}} and \gene{{PRRC2B}}. Variants altering the \gene{{alt-PRRC2A/B}} protein sequence ($n=23$) are highlighted in bold.\par}}
\vspace{{3pt}}
{{\scriptsize\setlength{{\tabcolsep}}{{3.2pt}}\renewcommand{{\arraystretch}}{{1.06}}
\begin{{longtable}}{{@{{}}llllllr@{{}}}}
\toprule
\textbf{{Gene}} & \textbf{{ClinVar ID}} & \textbf{{HGVS c.}} & \textbf{{ClinVar Significance}} & \textbf{{Frame-0 Effect}} & \textbf{{Frame-+2 Effect}} & \textbf{{+2 Pos}} \\
\midrule
\endfirsthead
\toprule
\textbf{{Gene}} & \textbf{{ClinVar ID}} & \textbf{{HGVS c.}} & \textbf{{ClinVar Significance}} & \textbf{{Frame-0 Effect}} & \textbf{{Frame-+2 Effect}} & \textbf{{+2 Pos}} \\
\midrule
\endhead
\bottomrule
\endfoot
{t4_lines}
\end{{longtable}}
}}

\newpage
\section*{{Supplementary Table S5: All 42 Single-Nucleotide Silent-Frame-0 / Nonsense-+2 Switches}}
\noindent{{\footnotesize Complete enumeration of every single-nucleotide substitution across human \gene{{PRRC2A}} ($n=22$) and \gene{{PRRC2B}} ($n=20$) that is strictly synonymous in Frame 0 while introducing a premature stop codon in Frame +2, audited directly against hg38 genomic DNA (including intronic flanks).\par}}
\vspace{{3pt}}
{{\fontsize{{6.2}}{{7.6}}\selectfont\setlength{{\tabcolsep}}{{2.2pt}}\renewcommand{{\arraystretch}}{{1.08}}
\begin{{longtable}}{{@{{}}lllllllll@{{}}}}
\toprule
\textbf{{Gene}} & \textbf{{CDS Sub.}} & \textbf{{hg38 Pos.}} & \textbf{{Ex.}} & \textbf{{Frame-0 (usage)}} & \textbf{{+2}} & \textbf{{PAM (p4--8)}} & \textbf{{Bystanders}} & \textbf{{Editor}} \\
\midrule
\endfirsthead
\toprule
\textbf{{Gene}} & \textbf{{CDS Sub.}} & \textbf{{hg38 Pos.}} & \textbf{{Ex.}} & \textbf{{Frame-0 (usage)}} & \textbf{{+2}} & \textbf{{PAM (p4--8)}} & \textbf{{Bystanders}} & \textbf{{Editor}} \\
\midrule
\endhead
\bottomrule
\endfoot
{t5_lines}
\end{{longtable}}
}}

\end{{document}}
"""
    out_tex = PREPRINT / "prrc2_plus2_supplement.tex"
    out_tex.write_text(tex)
    print(f"Wrote {out_tex} (35 screen hits, 10 lineages, 45 ClinVar, 42 switches)")


if __name__ == "__main__":
    build_supplement_tex()
