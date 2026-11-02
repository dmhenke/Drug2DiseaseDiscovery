# SLE drug discovery plan

Team notes from Yannis and David (received 2026-11-01), checked against what the
pipeline has already produced and against public sources. Every claim carries a
tag: **Known** (checked here, with a source), **Inferred** (follows from stated
assumptions), or **Guess** (unchecked).

This file sits outside `docs/` on purpose. If GitHub Pages is enabled, everything
under `docs/` becomes public, including Markdown files.

## Goal

Find small molecules that move the SLE transcriptional state of immune cells back
toward healthy.

**Focus requested by Yannis:** SLE patients who also have asthma (autoimmune plus
respiratory disease). See [Asthma subgroup](#asthma-subgroup) for whether the data
can support this.

**Method under consideration:** TranSiGen (Tong et al. 2024, *Nat Commun*,
[doi:10.1038/s41467-024-49620-3](https://www.nature.com/articles/s41467-024-49620-3),
[PMC11199551](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11199551/)). It is a deep
generative model. Given a cell line's basal expression and a compound's structure,
it predicts that compound's transcriptional profile, so it can score compounds that
were never assayed (**Known**, from the abstract). It is presumably trained on LINCS
L1000 perturbation data (**Inferred**; the full text has not been checked).

Two consequences for the team's questions (**Inferred**):
- **Candidate compounds need a structure (SMILES).** Antibodies (belimumab,
  anifrolumab, obinutuzumab, rituximab) cannot be scored and cannot serve as
  positive controls.
- **Candidate cell lines need a basal expression profile** available to the model.
  That limits the choice to lines already profiled in its training data.

## Where the project stands (2026-11-01)

| Step | Status | Where |
|---|---|---|
| SLE vs healthy signatures per cell type: Perez 2022 cohort 4, 20 v 20 donors, UCE embedding | Done, diagnostic. 12/12 interferon genes in classical monocytes | `docs/lupus-uce-test.html` |
| Drug matching against L1000 (L1000CDS2, reverse mode) | Done, exploratory. Only the monocyte signatures beat random gene sets | `docs/drug-reversal.html`, `scripts/l1000_reversal.py` |
| Full cohort-4 DE (all donors, all cells) | Not run. Recommended in PCA mode with author labels (about 3 days if run with UCE) | — |
| Second SLE study, cross-study design | Not started. Nehar-Belaid 2020 (GSE135779) is the top candidate | — |
| Bulk whole-blood meta-analysis (4 GEO microarray studies, limma), from the team's `sle_analysis` repo | Imported. Reviewed: the gene-intersection bug drops most interferon genes. Concordance with scRNA computed | `analysis/sle_bulk_meta/` |
| Chemical database and perturbation model design | ChEMBL chosen; seed table of 129 compounds built; chemCPA-style model on OP3 proposed | `planning/chemical_perturbation_model.md`, `shared/compounds/seed_compounds.tsv` |
| TranSiGen | Not started (bulk L1000 model; the single-cell design above is chemCPA-style instead) | — |

## Notes → actions

### Yannis

| # | Question | What we know now | Next step |
|---|---|---|---|
| Y1 | Which cell line to model? Data-driven approach? | The SLE signal is strongest in classical monocytes (68 DE genes; **Known**). Of the 100 monocyte hits in the L1000 lookup, only 2 came from myeloid lines (THP-1, SKM-1); most came from HA1E, VCaP, MCF7 and A549 (**Known**). L1000 is dominated by epithelial cancer lines. | Data-driven choice: rank candidate lines by how closely their basal expression matches the classical-monocyte pseudobulk profile. Restrict to lines that have a basal profile available to TranSiGen. THP-1 and U937 are the obvious candidates to test against that ranking (**Inferred**). |
| Y2 | Which drug database: ~30k CMap compounds, 600k natural products, or FDA-approved (needs SMILES)? | **Decided: ChEMBL 37 as master, Repurposing Hub as the curated subset, joined on InChIKey** (`planning/chemical_perturbation_model.md`). Seed table: 129 compounds, all with valid SMILES. | Start with FDA-approved plus CMap compounds. They have known mechanisms (DrugMechDB, ChEMBL), so hits can be interpreted. The 600k natural products come later: they are mostly out of the model's training domain, so predictions there are least reliable (**Inferred**). |
| Y3 | Other diseases to apply this to? Positive controls? | The paper validated TranSiGen in pancreatic cancer (**Known**, abstract). | Positive-control **compounds**: small-molecule SLE drugs (table below), whose predicted profiles should reverse the SLE signature better than random compounds do. Positive-control **diseases**: other interferon-driven conditions (dermatomyositis, Sjögren's, COVID-19). The same signature-reversal logic should recover interferon-pathway inhibitors there too (**Inferred**). |
| Y4 | List current SLE therapeutics | See [Current SLE therapeutics](#current-sle-therapeutics). | Get SMILES for the small molecules from ChEMBL. |

### David

| # | Item | What we know now | Next step |
|---|---|---|---|
| D1 | Deeper analysis of the lupus profile: GSEA against the literature; correlate with proteomics | **Bulk meta-analysis imported** (`analysis/sle_bulk_meta/`). Where single-cell hits are on the bulk table, they mostly replicate: 17/20 classical monocyte, 10/12 non-classical, 8/9 CD4 T, 4/12 NK (**Known**). The bulk table is missing most interferon genes because of a gene-intersection bug (README, problem 1). 26/31 ribosomal genes are also down in bulk blood. The GSEA draft doesn't run yet. The proteomics dataset listed (D5.12) contains no SLE cells (see below). | Run pre-ranked GSEA (fgsea, Hallmark and Reactome) per cell type. Compare with the interferon signature in Perez 2022 and Nehar-Belaid 2020. For protein-level data, use CITE-seq of SLE patients. None is in hand yet. |
| D2 | Show single genes or pathways in the top drug predictions. Is IFI27 down-regulated? | **IFI27 is up in SLE**, not down: log2FC 3.74 in classical monocytes (padj 7e-5), 5.20 in non-classical monocytes, 3.66 in CD8 T. It is also the top gene in the bulk meta-analysis (log2FC 4.59, adj.P 6e-79) (**Known**). Among the top 50 L1000 matches, IFI27 is pushed down by the compound in 20/50 (classical) and 25/50 (non-classical) signatures (**Known**). | The drug page now has a gene filter: type a gene to see which compounds reverse it. Pathway-level views come after GSEA (D1). |
| D3 | Find scRNA-seq data from SLE patients | See [Dataset registry](#dataset-registry). | Next study: GSE135779. |
| D4 | Drug/chemical databases: DrugMechDB, ChEMBL, Drug Hunter list | DrugMechDB appears twice in the notes; it's one resource. | ChEMBL for structures and targets. DrugMechDB for mechanism paths that explain a hit. LINCS L1000 is the perturbation reference already in use. |
| D5 | Cell lines: THP-1 (monocytic), a lymphocyte line (B, or preferably T), murine models | See Y1. No public mouse lupus spleen/blood scRNA-seq set with healthy controls was found (2026-11-01 search; **Known** for that search, not exhaustive). | T-cell line: Jurkat is the usual choice. Whether it has L1000 or TranSiGen coverage has not been checked (**Guess**). |
| D-note | TNFSF13B and OAS1 have the best diagnostic performance and are mainly expressed by monocytes ([PMC9627509](https://pmc.ncbi.nlm.nih.gov/articles/PMC9627509/)) | Our data agree. **TNFSF13B (BAFF)**: up in classical monocytes (log2FC 0.38, padj 0.027) and non-classical monocytes (0.74, padj 0.038). **OAS1**: up in both (0.98 and 0.89, padj < 0.02) (**Known**). BAFF is the target of belimumab. **None of the top 50 L1000 matches reverses TNFSF13B** (**Known**); whether L1000 measures it directly has not been checked. | Treat TNFSF13B as a second signature axis beside the interferon genes. A compound that reverses both is more interesting than one that only damps interferon genes. |

## Dataset registry

The notes list several entries that are **the same dataset**. Pooling them would
count donors two or three times.

| Notes entry | What it is | Verdict |
|---|---|---|
| D5.10 GSE174188 · D5.8 HCA "SLE Immune Profiling" · D6.1 Perez 2022 *Science* | Perez et al. 2022 (PMID 35389781): 1.26M PBMCs, 162 SLE / 99 healthy. Also on CELLxGENE (dataset 218acb0f). | **Already in use.** One dataset under three names (**Known**: same counts and paper). |
| D5.7 / D5.11 GSE137029 (PRJNA564301) | Earlier release from the same paper: ~1M PBMCs, 134 SLE / 58 healthy, plus 15 flare cases. The paper's data statement cites GSE137029. | **Overlaps Perez heavily**: same paper and cohorts (CLUES, ImmVar) (**Inferred**; donor-level overlap not checked). Do not add it as a second study. |
| D5.5 GSE162577 | SLE PBMC from China. The notes say 2 SLE + 1 healthy; a 2024 review lists 1 SLE + 1 control. | **Too small for donor-level DE** (n ≤ 2 per group). Skip. |
| D5.4 GSE163121 | "B-cell profiling in SLE" | **Unverified**: the accession didn't come up in search. Check GEO before planning around it. |
| D5.6 GSE186476 | Lupus skin, lesional and non-lesional, plus blood; scRNA-seq and spatial (PMID 35290245) | Real. **Different tissue.** Useful for a skin-specific question, not for the blood signature. The notes also mention kidney; not confirmed. |
| D5.9 SLE-ImmuneAtlas (gao-lab, GitHub) | scRNA and scVDJ | **Unverified.** Find which GEO studies it combines before using it; they may include Perez. |
| D5.12 scProteomics "Broad immune activation…" | Kotliarov et al. 2020 *Nat Med* (doi:10.1038/s41591-020-0769-8). CITE-seq of 53,201 cells from **healthy influenza-vaccine responders**. The lupus link comes from separate SLE data. | **Contains no SLE cells** (**Known**). Not usable as SLE proteomics. |
| Not in the notes | Nehar-Belaid et al. 2020 (GSE135779): childhood SLE PBMC, 33 SLE / 11 healthy, plus an adult validation set | **Recommended next study.** Independent cohort with both conditions; tests the `~ study + condition` design. |

## Current SLE therapeutics

These are the positive-control candidates. Only small molecules can be scored by
TranSiGen or L1000.

| Drug | Target / mechanism | US status for SLE (checked 2026-11-01) | Small molecule? |
|---|---|---|---|
| Hydroxychloroquine | Antimalarial; TLR7/9 and endosomal signaling | Standard of care | Yes |
| Glucocorticoids (prednisone, methylprednisolone) | Glucocorticoid receptor | Standard of care | Yes |
| Mycophenolate mofetil | IMPDH inhibitor | Used off-label, standard in lupus nephritis | Yes |
| Azathioprine, methotrexate, cyclophosphamide | Antimetabolite / alkylating | Used | Yes |
| Voclosporin | Calcineurin inhibitor | Approved for lupus nephritis, Jan 2021 (**Known**) | Cyclic peptide, outside the usual chemical space |
| Belimumab | Anti-BAFF (TNFSF13B) | Approved, SLE and lupus nephritis | No (antibody) |
| Anifrolumab | Anti-IFNAR1 | Approved for SLE. Subcutaneous form approved April 2026 (**Known**) | No |
| Obinutuzumab | Anti-CD20 | Approved for lupus nephritis, Oct 2025. SLE application accepted April 2026, decision expected by Dec 2026 (**Known**) | No |
| Rituximab | Anti-CD20 | Off-label | No |
| Deucravacitinib | TYK2 inhibitor | Phase 3 in SLE ongoing; no phase 3 readout found (**Known**, 2026-11-01) | Yes |
| JAK inhibitors (baricitinib, tofacitinib) | JAK1/2, JAK1/3 | Not approved for SLE; baricitinib's phase 3 program did not succeed (**Guess**, from memory) | Yes |

**Positive-control test (proposed):**
- Hydroxychloroquine, dexamethasone or prednisolone, mycophenolic acid and
  deucravacitinib should score better at reversing the monocyte signature than
  random compounds from the same library.
- An interferon-pathway inhibitor (TYK2 or JAK) should specifically reverse the
  interferon genes.
- If none of these beat random, the drug-matching method is not informative for SLE.

## Asthma subgroup

- **Epidemiology:**
  - Asthma is associated with higher SLE risk: pooled OR 1.37 (95% CI 1.14–1.65), from 10 studies with high heterogeneity (Charoenngam et al. 2020, *Clin Rheumatol*).
  - In the FORWARD SLE cohort, 19.8% reported asthma at entry, and asthma tracked with worse patient-reported outcomes (Katz et al. 2021, *ACR Open Rheumatol*) (**Known**).
- **Data gap:** none of the SLE single-cell datasets above records asthma status.
  - Perez 2022's metadata has no comorbidity fields (**Known**, from the obs columns).
  - The others have not been checked (**Guess**: unlikely).
- **Feasible version:** find genes and pathways shared between the SLE blood signature and an asthma blood or airway signature, from separate datasets.
  - Rank compounds that reverse both.
  - Flag compounds known to worsen asthma (e.g. non-selective beta-blockers).
- **What it can't do:** identify mechanisms specific to patients who have both diseases. That needs a cohort with both measured.
- **Lung data:** the Human Lung Cell Atlas on CELLxGENE covers lung disease but not asthma among its listed diseases (**Known**, CELLxGENE listing). An asthma scRNA-seq source still needs to be found.

## Open questions for the team

1. Is the target population SLE generally, or specifically SLE with asthma? That decides whether the next dataset is a second SLE study (GSE135779) or an asthma dataset.
2. Small molecules only? If so, the approved antibody drugs drop out as positive controls, and the control list is the small-molecule table above.
3. Cell line for TranSiGen: accept a data-driven ranking (Y1), or fix THP-1 in advance?
4. Should TranSiGen run on the node? It has an 8 GB GPU. Model size and requirements are not yet checked.

## Sources

- Tong et al. 2024, TranSiGen: [Nat Commun](https://www.nature.com/articles/s41467-024-49620-3) · [PMC11199551](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11199551/)
- Perez et al. 2022, *Science* 376:eabf1970, PMID 35389781 · [GSE174188](https://www.omicsdi.org/dataset/geo/GSE174188) · [GSE137029](https://omicsdi.org/dataset/geo/GSE137029)
- [GSE186476](https://www.omicsdi.org/dataset/geo/GSE186476) · [GSE162577 in a 2024 review table](https://pmc.ncbi.nlm.nih.gov/articles/PMC11751531/table/tbl5)
- Kotliarov et al. 2020, [Nat Med](https://www.nature.com/articles/s41591-020-0769-8.pdf)
- [Obinutuzumab approved for lupus nephritis](https://www.the-rheumatologist.org/article/fda-approves-obinutuzumab-for-lupus-nephritis/) · [SLE application accepted](https://www.roche.com/investors/updates/inv-update-2026-04-21)
- [Voclosporin approval](https://hcplive.com/view/fda-approves-voclosporin-lupus-nephritis) · [Deucravacitinib phase 3, NCT05620407](https://clinicaltrials.gov/study/NCT05620407)
- [Lupus Foundation: treatments being studied](https://www.lupus.org/resources/treatments-being-studied-for-lupus)
- Asthma and SLE: [Katz et al. 2021](https://doaj.org/article/7eb3c4d6c6764e2993878d1e70f202f0) · [Charoenngam et al. 2020 meta-analysis](https://repository.li.mahidol.ac.th/handle/123456789/58281)
- TNFSF13B and OAS1: [PMC9627509](https://pmc.ncbi.nlm.nih.gov/articles/PMC9627509/)
- DrugMechDB: [Sci Data 2023](https://www.nature.com/articles/s41597-023-02534-z) · [ChEMBL](https://www.ebi.ac.uk/chembl/) · [Drug Hunter database list](https://drughunter.com/resource/drug-discovery-websites-and-databases)
