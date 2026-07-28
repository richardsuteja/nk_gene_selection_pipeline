# Literature Validation — amiRNA Knockdown Candidates

Systematic search for direct causal in vivo NK or CD8+ T cell knockout/knockdown
evidence per candidate, run against the 19-TF list from `06_centrality_VIPER.R`
(`S06_results18_amiRNA_knockdown_candidates_VIPER.csv`).

**Result Selection Criteria:** 
1. Indexed in PubMed
2. Single-target KO/KD in NK/CD8+ T cells
3. Available results from in vivo assay

## Summary

| Gene | Search summary | Verdict |
|---|---|---|
| **NR4A1** | 29/07/2026, 114 studies. 3 in vivo NK/CD8+ findings, all positive (Yu 2023; Nowyhed 2015 ×2). | **Included** |
| **ATF3** | 29/07/2026, 7 studies. 1 in vivo NK finding, positive (Rosenberger 2008). | **Included** |
| FOSB | 29/07/2026, 45 studies. No qualifying in vivo evidence. | Lack of information |
| FOS | 29/07/2026, 68 studies. No qualifying in vivo evidence. | Lack of information |
| ZNF683 | 29/07/2026, 11 studies. No qualifying in vivo evidence. | Lack of information |
| YBX1 | 29/07/2026, 10 studies. No qualifying in vivo evidence. (A CD8+ T cell in vivo shRNA-KD finding exists in the literature but is not PubMed-indexed — fails the inclusion criterion outright.) | Lack of information |
| MXD1 | 29/07/2026, 1 study. No qualifying in vivo evidence. | Lack of information |
| REL | 29/07/2026, 17 studies. 2 neutral in vivo findings (Tato 2007; Liou 1999); Saibil 2007 reports a 3× in vitro viability decrease. | Excluded — in vitro negative tendency |
| EGR3 | 29/07/2026, 37 studies. No qualifying in vivo evidence. | Lack of information |
| GPBP1 | 29/07/2026, 0 studies. No qualifying in vivo evidence. | Lack of information |
| ZFHX2 | 29/07/2026, 0 studies. No qualifying in vivo evidence. (Not to be confused with ZHX2, a distinct gene.) | Lack of information |
| JUNB | 29/07/2026, 52 studies. 1 negative finding (Sarkar 2025). | Excluded |
| JUN | 29/07/2026, 91 studies. No qualifying in vivo evidence. | Lack of information |
| EGR1 | 29/07/2026, 18 studies. 1 negative finding (Singh 2004). | Excluded |
| NR4A2 | 29/07/2026, 410 studies. No qualifying in vivo evidence. | Lack of information |
| ZNF331 | 29/07/2026, 1 study. No qualifying in vivo evidence. | Lack of information |
| TOPORS | 29/07/2026, 0 studies. No qualifying in vivo evidence. | Lack of information |
| KLF6 | 29/07/2026, 5 studies. No qualifying in vivo evidence. | Lack of information |
| MYB | 29/07/2026, 22 studies. 1 negative finding (Gautam 2019) — reduced survival, memory formation and self-renewal despite increased GZMB on primary stimulation. | Excluded |

**Result: NR4A1 and ATF3 are the only two candidates with direct causal in vivo**

## Detailed findings (genes with specific in vivo results)

| Gene | Author (PMID) | Outcome | Verdict |
|---|---|---|---|
| NR4A1 | Yu (2023), 36420610 | CRISPR-mediated KO in NK cells upregulates cytotoxic markers (CD226, CD107a, granzyme B, IFN-γ) and increases anti-PD1 responsiveness | Positive |
| NR4A1 | Nowyhed (2015), 25762306 | Deficiency/KD in CD8+ T cells upregulates Runx3 and its target genes (IFNG, EOMES, GZMB, PRF1); increases mature CD8+ T cell population | Positive |
| NR4A1 | Nowyhed (2015), 26363057 | Deficiency in CD8+ T cells upregulates IRF4; increases post-stimulation proliferation | Positive |
| ATF3 | Rosenberger (2008), 18268321 | Deficiency in NK cells increases IFN-γ mRNA; no change to perforin/granzyme B/TNF/GM-CSF mRNA | Positive |
| REL | Tato (2007), 16481345 | Deficiency in NK cells does not alter post-stimulation IFN-γ or proliferation | Neutral |
| REL | Liou (1999), 10221648 | Deficiency in CD8+ T cells does not impair thymic distribution | Neutral |
| REL | Saibil (2007) | In vitro: 3× viability decrease (24.8%→71.7%) in CD8+ T cells vs. WT | Negative (in vitro) |
| JUNB | Sarkar (2025), 39425978 | Deficiency in CD8+ T cells inhibits post-stimulation clonal expansion; reduces ID2/TBX21 | Negative |
| EGR1 | Singh (2004), 15356133 | Deficiency in CD8+ T cells does not alter IFN-γ/TNF-α/IL-2R; reduces post-stimulation clonal expansion | Negative |
| MYB | Gautam (2019), 30778251 | Deficiency in CD8+ T cells lowers survival, memory formation and self-renewal (stemness), reduces oxidative metabolism and splenic accumulation — despite increased GZMB on primary stimulation | Negative |

## Search strings

Timestamp 29/07/2026 for all.

**FOSB**
```
("FOSB" OR “FOS-B” OR “G0S3” OR “GOSB” OR “GOS3” OR “AP-1” OR “MGC42291” OR “DKFZp686C0818" OR "FBJ Murine Osteosarcoma Viral Oncogene Homolog B” OR “G0/G1 Switch Regulatory Protein 3”) AND ("Knockout" OR "KO" OR "Knockdown" OR "KD") AND ("CD8" OR "CD8+" OR "NK" OR "Natural Killer")
```

**FOS**
```
("FOS" OR "c-FOS" OR "AP-1" OR "FBJ murine osteosarcoma viral oncogene homolog" OR "Activator Protein 1") AND ("Knockout" OR "KO" OR "Knockdown" OR "KD") AND ("CD8" OR "CD8+" OR "NK" OR "Natural Killer")
```

**ZNF683**
```
("ZNF683" OR "HOBIT" OR “MGC33414” OR “Zinc Finger Protein 683” OR “Homolog of Blimp-1”) AND ("Knockout" OR "KO" OR "Knockdown" OR "KD") AND ("CD8" OR “CD8+” OR "NK" OR “Natural Killer”)
```

**YBX1**
```
("YBX1” OR “YB-1” OR “YB1” OR “DBPB” OR “MDR-NF1” OR “BP-8” OR “CSDB” OR “CSDA2” OR “NSEP1” OR “Y-box binding protein 1” OR “DNA-binding protein B” OR “Multidrug resistance-associated protein / neurofibromatosis 1 associated“ OR “Binding protein 8” OR “Cold shock domain-containing protein A2” OR “Cold shock domain-containing protein B” OR “nuclease sensitive element binding protein 1”) AND ("Knockout" OR "KO" OR "Knockdown" OR "KD") AND ("CD8" OR “CD8+” OR "NK" OR “Natural Killer”)
```

**MXD1**
```
("MXD1" OR “bHLHc58” OR “MAX dimerization protein 1”) AND ("Knockout" OR "KO" OR "Knockdown" OR "KD") AND ("CD8" OR "CD8+" OR "NK" OR “Natural Killer”)
```

**REL**
```
("REL” OR “I-Rel” OR “c-REL” OR “HIVEN86A”) AND ("Knockout" OR "KO" OR "Knockdown" OR "KD") AND ("CD8" OR “CD8+” OR "NK" OR “Natural Killer”)
```

**EGR3**
```
("EGR3” OR “PILOT” OR “Early Growth Response 3”) AND ("Knockout" OR "KO" OR "Knockdown" OR "KD") AND ("CD8" OR "CD8+" OR "NK" OR “Natural Killer”)
```

**GPBP1**
```
("GPBP1” OR “GC-rich promoter binding protein 1” OR “Vasculin” OR “DKFZp761C169”) AND ("Knockout" OR "KO" OR "Knockdown" OR "KD") AND ("CD8" OR "CD8+" OR "NK" OR “Natural Killer”)
```

**ZFHX2**
```
("ZFHX2” OR “ZNF409” OR “KIAA1762” OR “KIAA1056” OR “ZFH-5” OR “zinc finger homeobox 2” OR “Zinc Finger Protein 409”) AND ("Knockout" OR "KO" OR "Knockdown" OR "KD") AND ("CD8" OR "CD8+" OR "NK" OR “Natural Killer”)
```

**JUNB**
```
("JUNB” OR “JUN-B” OR “AP-1”) AND ("Knockout" OR "KO" OR "Knockdown" OR "KD") AND ("CD8" OR "CD8+" OR "NK" OR “Natural Killer”)
```

**JUN**
```
("c-JUN” OR “v-JUN Avian Sarcoma Virus 17 Oncogene Homolog” OR “AP-1”) AND ("Knockout" OR "KO" OR "Knockdown" OR "KD") AND ("CD8" OR "CD8+" OR "NK" OR “Natural Killer”)
```

**ATF3**
```
("ATF3” OR “Activating Transcription Factor 3”) AND ("Knockout" OR "KO" OR "Knockdown" OR "KD") AND ("CD8" OR “CD8+” OR "NK" OR “Natural Killer”)
```

**EGR1**
```
("EGR1” OR “TIS8” OR “AT225” OR “G0S30” OR “NGFI-A” OR “ZIF268” OR “KROX-24” OR “ZIF-268” OR “Early Growth Response 1” OR “Early Growth Response Protein 1” OR “TPA-Induced Sequence 8” OR “Zinc Finger Gene 225” OR “G0/G1 Switch Regulatory Gene 30” OR “Nerve Growth Factor-Induced Protein A” OR “Transcription Factor ETR103”) AND ("Knockout" OR "KO" OR "Knockdown" OR "KD") AND ("CD8" OR "CD8+" OR "NK" OR “Natural Killer”)
```

**NR4A2**
```
("NR4A2” OR “NURR1” OR “TINUR” OR “RNR1” OR “HZF3” OR “Nuclear Receptor Subfamily 4 Group A Member 2” OR “Nuclear Receptor-Related 1” OR “Transcriptionally-Inducible Nuclear Receptor” OR “Nuclear receptor of T cells” OR “Nuclear Receptor-Related 1” OR “Heart Zinc Finger 3”) AND ("Knockout" OR "KO" OR "Knockdown" OR "KD") AND ("CD8" OR “CD8+” OR "NK" OR “Natural Killer”)
```

**ZNF331**
```
("ZNF331” OR “ZNF463” OR “ZNF361” OR “RITA” OR “Zinc Finger Protein 331” OR “Zinc Finger Protein 463” OR “Zinc Finger Protein 361” OR “Rearranged in Thyroid Adenoma”) AND ("Knockout" OR "KO" OR "Knockdown" OR "KD") AND ("CD8" OR "CD8+" OR "NK" OR “Natural Killer”)
```

**TOPORS**
```
("TOPORS” OR “RP31” OR “TP53BPL” OR “TOP1 binding arginine/serine rich protein” OR “retinitis pigmentosa 31” OR “Topoisomerase I binding arginine/serine-rich E3 ubiquitin protein ligase”) AND ("Knockout" OR "KO" OR "Knockdown" OR "KD") AND ("CD8" OR "CD8+" OR "NK" OR “Natural Killer”)
```

**KLF6**
```
("KLF6” OR “BCD1” OR “ST12” OR “COPEB” OR “CPBP” OR “GBF” OR “Zf9” OR “PAC1” OR “Kruppel Like Factor 6” OR “B-Cell-Derived protein 1” OR “Core Promoter Element Binding Protein” OR “Suppressor of Tumorigenicity 12” OR “GC-Rich Sites-Binding Factor”) AND ("Knockout" OR "KO" OR "Knockdown" OR "KD") AND ("CD8" OR “CD8+” OR "NK" OR “Natural Killer”)
```

**NR4A1**
```
("NR4A1” OR “HMR” OR “GRFP1” OR “TR3” OR “N10” OR “NAK-1” OR “NGFIB” OR “NUR77” OR “Nuclear Receptor Subfamily 4 Group A Member 1” OR “Glucocorticoid Receptor-Interacting Protein 1” OR “Testicular Receptor 3” OR “Nerve Growth Factor-Induced Gene B” OR “Nuclear Orphan Receptor 1” OR “Growth Factor-Inducible Nuclear Protein N10”) AND ("Knockout" OR "KO" OR "Knockdown" OR "KD") AND ("CD8" OR “CD8+” OR "NK" OR “Natural Killer”)
```

**MYB**
```
("MYB” OR “c-MYB” OR “v-myb avian myeloblastosis viral oncogene homolog”) AND ("Knockout" OR "KO" OR "Knockdown" OR "KD") AND ("CD8" OR “CD8+” OR "NK" OR “Natural Killer”)
```
