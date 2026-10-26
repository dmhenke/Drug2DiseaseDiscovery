# Drug2DiseaseDiscovery

Cross-study, cross-species single-cell RNA-seq integration and differential
expression built on **UCE** (Universal Cell Embedding; Rosen *et al.*, *Nature*
2026). UCE maps cells from any study/platform/species into one shared space so
they sit on the same footing; differential expression then runs on expression
counts with `study` as a batch covariate (DESeq2 `~ study + condition`).

**Lung is the first study.** The project is organized so that adding other
organs / illnesses later is a copy-one-folder operation that never touches or
overwrites existing work.

## Folder layout

```
Drug2DiseaseDiscovery/
├── README.md                     ← you are here
├── docs/
│   └── UCE_implementation_guide.md   the method, rationale, and run order
├── pipeline/                     ← reusable, ORGAN-AGNOSTIC code (edit rarely)
│   ├── setup_uce_gpu.sh              one-time GPU/UCE install (shared)
│   ├── 01_build_species_anndata.R    Stage 0: assemble per-species h5ads
│   ├── run_uce.sh                    Stage 1: UCE embedding (Python/GPU)
│   ├── 02_integrate_and_cluster.R    Stage 2: joint clustering + annotation
│   ├── 03_pseudobulk_de.R            Stage 3: pseudobulk DESeq2
│   └── UCE/                          (created by setup; the cloned UCE repo)
├── shared/                       ← reused by every study
│   ├── models/                       UCE weights (download once)
│   └── references/                   reference atlases for label transfer
└── studies/                      ← one self-contained folder PER study
    ├── _template/                    copy this to start a new organ/illness
    │   ├── config.R                  study-specific settings
    │   ├── data/                     raw input .h5ad files
    │   ├── interim/                  intermediate + UCE-embedded files
    │   └── results/                  DE tables, warnings
    └── lung/                         the first study
        ├── config.R
        ├── data/    interim/    results/
```

The rule: **code is shared and generic; data, config, and outputs are
per-study.** Because each study writes only inside its own `interim/` and
`results/`, running a new organ can never overwrite another.

## Running a study (lung)

Prereqs once: `pipeline/setup_uce_gpu.sh` (GPU box), and the R packages listed
in `docs/UCE_implementation_guide.md` §5. Put raw lung `.h5ad` files in
`studies/lung/data/` and list them in `studies/lung/config.R`. Then:

```bash
Rscript pipeline/01_build_species_anndata.R studies/lung
pipeline/run_uce.sh                          studies/lung     # GPU step
Rscript pipeline/02_integrate_and_cluster.R  studies/lung
Rscript pipeline/03_pseudobulk_de.R          studies/lung
```

Results appear in `studies/lung/results/` (per-cell-type DE tables, a combined
table per species, a human–mouse conserved-ortholog table, and `warnings.txt`).

## Adding another organ or illness

```bash
cp -r studies/_template studies/kidney      # or liver_nafld, retina_amd, ...
# edit studies/kidney/config.R  (STUDY_NAME, INPUTS, CONDITION_MAP, ...)
# drop raw .h5ad files into studies/kidney/data/
Rscript pipeline/01_build_species_anndata.R studies/kidney
pipeline/run_uce.sh                          studies/kidney
Rscript pipeline/02_integrate_and_cluster.R  studies/kidney
Rscript pipeline/03_pseudobulk_de.R          studies/kidney
```

No pipeline code changes; the shared UCE model and any reference atlases in
`shared/` are reused. Each study is independent.

## Notes

- UCE embedding (Stage 1) is the only Python/GPU step; everything else is R.
- Use the 33-layer model throughout — its embeddings are not interchangeable
  with the 4-layer model, so keep one study's embeddings on a single model.
- Sources: paper https://www.nature.com/articles/s41586-026-10689-z ·
  code https://github.com/snap-stanford/UCE
