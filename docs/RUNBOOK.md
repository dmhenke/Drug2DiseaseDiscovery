# Runbook: test on this Mac, deploy the full cohort on a GPU

Status tags used below:
- **Tested** = run here (Mac, M4, 16 GB, 2026-10-25).
- **Untested** = written, not executed (needs R packages, UCE, or a GPU).
- **Inferred** = from reading the code or from arithmetic. **Guess** = a number I could not measure.

## 0. Approvals needed before any of this

These steps download or install things from outside the machine. Each needs your yes first.

| Step | What | Size / source |
|---|---|---|
| R packages | Bioconductor and CRAN packages, incl. `RcppAnnoy` | Bioconductor, CRAN |
| UCE code | `git clone` of `snap-stanford/UCE` (`scripts/setup_uce_mac.sh`) | GitHub |
| Python env | `pip` installs for UCE's pins | PyPI |
| UCE weights | `33l_8ep_1024t_1280.torch` into `shared/models/` | figshare. Size unknown: the HEAD request returned no length |
| UCE support files | Token file, ESM2 embeddings, offsets. UCE downloads these **automatically on its first run** (`evaluate.py`, `check_paths`) | figshare; sizes unknown |
| Test data | CELLxGENE trial collection, T Cells file | 1.08 GB, `datasets.cellxgene.cziscience.com` |
| Full data | Perez et al. PBMC object | 12.2 GB, `datasets.cellxgene.cziscience.com` |

## 1. Hardware and time expectations

| Item | This Mac (M4, 16 GB) | GPU host (recommended) |
|---|---|---|
| Accelerator | Metal (MPS), **untested** with UCE | NVIDIA A100 or H100, 80 GB |
| Host RAM | 16 GB: tight for the 33-layer model load | ≥ 128 GB (Stage 0 and 2 hold the full cohort in memory) |
| Disk free | ≥ 60 GB for a 20k-cell test | ≥ 400 GB for the full cohort |
| UCE batch | 4 (**Inferred**: attention activations scale with batch and sequence length 1,536) | 25 (paper's setting) |

**Time (Inferred from code; throughput is a Guess):**
- Per cell, the 33-layer model costs about 2 TFLOP (**Inferred** from `evaluate.py`: 1,280-dim model, 5,120 hidden, 33 layers, sequence length 1,536).
- This Mac: about 1–2 TFLOP/s effective (**Guess**). That gives roughly **3–6 h per 10k cells** (**Inferred** from the first two).
- Full cohort (1.26M cells) on an A100: 20–40 TFLOP/s effective (**Guess**), so **about 18–35 h** for Stage 1.
- **Correction:** the 4–6 h GPU figure in earlier reports was wrong. It does not follow from the 1.5 h / 80k-cell estimate, which scales to about a day for 1.26M cells.
- The benchmark (Part B, step B2) replaces every Guess here with a measured rate.

The 4-layer model is about 7× cheaper but is a **different embedding**. The README says the two are not interchangeable. Use it only to test plumbing, and never to report results.

## Part A. Test on this Mac

### A1. Preflight (read-only, **tested**)

```bash
scripts/preflight.sh studies/synthetic_t0
```

Expected here now: R packages FAIL until A2 is done, and UCE FAIL until section 0 is approved. Python and disk should pass. `run_study.sh` aborts at preflight until the R packages are installed.

### A2. Install the R packages (**untested**, needs approval)

```bash
Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager"); BiocManager::install(c("zellkonverter","SingleCellExperiment","S4Vectors","scuttle","DESeq2","bluster","BiocNeighbors","RcppAnnoy","igraph","Matrix"))'
```

### A3. Synthetic fixture: no downloads, no UCE (**tested** to the point of writing files; R stages untested)

The fixture is already built in `studies/synthetic_t0/` (1,980 cells, 3 studies, 40 planted genes).
Rebuild it if you change the seed:

```bash
python3 scripts/make_synthetic_fixture.py
```

Run the full pipeline using PCA in place of UCE. This tests Stages 0, 2, and 3 with no model:

```bash
UCE_MODE=pca scripts/run_study.sh studies/synthetic_t0 diagnostic
```

Then check the DE results against the planted truth:

```bash
python3 scripts/check_t0.py studies/synthetic_t0
```

Expectations, written before the run: recall of planted genes ≥ 0.75 in every cell type, and ≤ 25 non-planted calls per cell type. The checker is **tested** on fabricated results only, so a PASS here proves nothing about the pipeline until the run itself passes.

Also run the T0 negative control (healthy vs healthy, see section 3). It is **not yet written**.

### A4. Real-data subsample (**tested** on the fixture; untested on real files)

Download the trial T Cells file only after approval (section 0). Inspect it first:

```bash
python3 scripts/subsample_h5ad.py data/cellxgene/trial_tcells.h5ad studies/lupus_test/data/lupus_test.h5ad --donor-col donor_id --condition-col disease --healthy normal --diseased "systemic lupus erythematosus" --donors-per-condition 5 --max-cells-per-donor 1000 --dry-run
```

The column names above are CELLxGENE conventions and **unverified** for this file. If the dry run errors, read `obs` with `anndata` and fix the names. Then run it for real (drop `--dry-run`). It writes raw counts to `X` and a manifest next to the output.

Edit `studies/lupus_test/config.R` to match the column names, then:

```bash
UCE_MODE=pca scripts/run_study.sh studies/lupus_test diagnostic
```

This plumbs real data through Stages 0–3 with no model.

### A5. UCE on the Mac (**untested**, needs approval)

```bash
scripts/setup_uce_mac.sh
```

Then point the scripts at the venv and check the device:

```bash
export UCE_PYTHON="$HOME/uce-venv/bin/python"
scripts/preflight.sh studies/lupus_test
```

Expected: MPS. If `preflight` reports `cpu`, stop. The 33-layer model is not practical on CPU. If the torch 2.1.1 pin fails on Apple Silicon, record the error and go to the GPU plan. Do not swap in an unpinned torch silently.

### A6. Benchmark before any real run (**untested**; the arithmetic is **tested**)

```bash
BATCH=4 scripts/benchmark_uce.sh studies/lupus_test/interim/lupus_test_human.h5ad 200
```

(Run Stage 0 first, or point at any species-level h5ad with raw counts.) The script prints the marginal rate and the estimated hours for 10k, 20k, and the full cohort. Use those numbers, not section 1, to plan.

### A7. Real UCE run on the test subsample (**untested**)

```bash
BATCH=4 SHARD_CELLS=5000 UCE_PYTHON="$HOME/uce-venv/bin/python" scripts/run_study.sh studies/lupus_test diagnostic
```

`SHARD_CELLS=5000` keeps the dense UCE preprocessing small (**Inferred** from memory use). Watch the log for the `RUN KIND` line at the top and bottom.

## Part B. Deploy the full cohort on a GPU

### B0. Provision

- One NVIDIA A100 or H100 with 80 GB, CUDA 12.1 driver (matches `setup_uce_gpu.sh`).
- Host RAM ≥ 128 GB. Disk ≥ 400 GB.
- Run `pipeline/setup_uce_gpu.sh` once (it installs conda, torch with CUDA, and UCE). It is the existing script and is **untested** by me.

### B1. Preflight

```bash
UCE_PYTHON=/path/to/uce/env/bin/python scripts/preflight.sh studies/lupus_full
```

Expect `accelerator: cuda`.

### B2. Benchmark on the GPU and record the decision

```bash
BATCH=25 scripts/benchmark_uce.sh studies/lupus_full/interim/lupus_full_human.h5ad 400
```

Stop here if the estimate for 1.26M cells is more than you can afford. Record the estimate in the run log.

### B3. Build the full input

Subsample with no donor cap, using every donor, and raw counts:

```bash
python3 scripts/subsample_h5ad.py data/cellxgene/perez_pbmc.h5ad studies/lupus_full/data/lupus_full.h5ad --donor-col donor_id --condition-col disease --healthy normal --diseased "systemic lupus erythematosus" --all-donors --max-cells-per-donor 0 --counts raw --dry-run
```

Check the dry run, then remove `--dry-run`. Confirm the donor counts and the `obs` names before the real run.

### B4. Full run

```bash
BATCH=25 SHARD_CELLS=50000 UCE_PYTHON=/path/to/uce/env/bin/python scripts/run_study.sh studies/lupus_full full
```

- The full run is slow (section 1). Run it inside `tmux` or `nohup` so a dropped SSH session does not kill it.
- `run_study.sh` does not resume. A failure restarts from Stage 0. Keep the log.
- Shard outputs are in `studies/lupus_full/interim/uce/`. The merge step fails loudly if any shard is missing or the cell order differs.

### B5. What to check after the run

- `results/run_full_*.log` has the RUN KIND line at the top and bottom, the code hash, and the parameters.
- `results/warnings.txt` lists skipped cell types and confounded designs.
- For lupus, check the interferon genes (`IFI44L`, `IFI6`, `ISG15`, `MX1`, `OAS1`, `IFIT1`, `IFIT3`, `LY6E`, `STAT1`) are up in SLE in most cell types. Compare with the report's T2a.
- Check `XIST` and `RPS4Y1` are not among the top hits (sex leakage, T2e). A design with no sex covariate will show them.

## 2. Data requirements

| | Synthetic T0 | Mac test (A4) | Full cohort (B3) |
|---|---|---|---|
| Format | h5ad, raw integer counts in `X` | same | same |
| Gene IDs | symbols, or Ensembl plus `feature_name` | same | same |
| Cells | 1,980 | ~10k | 1.26M |
| Donors | 18 (3 per condition per study) | 5 per condition | all |
| Studies | 3 | 1 | 1 (add a second cohort later) |
| Labels | `donor`, `condition` (names vary) | `donor_id`, `disease`, verified by dry run | same |
| Disk | 8 MB | ~0.1–0.5 GB | ~12–20 GB |

## 3. Known gaps (not done)

- **T0 negative control** (healthy vs healthy, expect ~0 hits): not written.
- **Stage 0 `make.unique` on duplicate symbols** (report C4): not fixed. Duplicate counts are discarded, not summed.
- **Cell-type confounding guard and results location** (C5): not fixed.
- **Shrunken fold changes** (S2): not added. Effect sizes in the output are still unshrunk MLEs.
- **Global multiplicity** (S1), **covariates** (S4), **donor-level aggregation** (S5), **rarefaction** (S6), **doublets and ambient RNA** (C7): not addressed.
- **Reproducibility across hardware:** seeds fix the sampling within one machine. MPS and CUDA can still differ numerically, so do not compare a Mac embedding with a GPU embedding.

## 4. Changes made to the pipeline

| File | Change | Report finding | Status |
|---|---|---|---|
| `pipeline/run_uce.sh` | Trailing slash on `--dir`; rename to the name Stage 2 reads; clear stale cache; seeded wrapper; PCA mode; shards | B1, B4, C2, T1 | Shell syntax checked. Stage 1 untested |
| `pipeline/01_build_species_anndata.R` | Sparse merge; intersect genes; unique cell IDs with an assert | B2, B3, C1 | Parse checked. Untested |
| `pipeline/02_integrate_and_cluster.R` | `set.seed`; Annoy k-NN | C2, S7 | Parse checked. Untested |
| `pipeline/03_pseudobulk_de.R` | Removed the dense `as.matrix` | B2 | Parse checked. Untested |

The R changes are untested. The `BNPARAM = BiocNeighbors::AnnoyParam()` argument name in `makeSNNGraph` should be checked against your installed `bluster` version.
