# Parallelized forecasting pipelines

Notes on the parallelization of the crime and macro forecasting exercises, the validation sequence used to check it, and the files involved.

------------------------------------------------------------------------

## 1. What changed

Both forecasting exercises loop over train/test splits, refitting every model at each split. Nothing carries state from one split to the next, so the loops are embarrassingly parallel. They now run through `future_lapply` instead of `for`.

At the original MCMC settings:

| Exercise            | Sequential          | 13 workers |
|---------------------|---------------------|------------|
| Crime, 289 origins  | \~90 min            | 10 min     |
| Macro, 116 vintages | \~2.9 hr (inferred) | 13.5 min   |

Scaling is linear in worker count, as expected. The practical gain is headroom rather than wall clock: the exercises can now be run at proper MCMC settings (5000 draws after 5000 burn-in rather than 1000 with none) in about an hour, where sequentially those settings would take most of a day. That is what made the factor-dimension sweeps feasible.

------------------------------------------------------------------------

## 2. How it works

The same pattern in both pipelines.

**The loop body became a function** of the split index that returns its own forecast slices as a list. It does not touch the large preallocated arrays, which matters because `future` would otherwise try to export them to every worker.

**The parent assembles.** After `future_lapply` returns, a short loop writes each result into the preallocated arrays by index. Scoring downstream is unchanged from the original scripts.

**Settings are passed explicitly** in a `settings` list rather than found by scoping. Workers do not see the calling environment's globals reliably.

**Workers source the helpers themselves** on first use, guarded by `exists("DGFC.mcmc")`. `future` does try to detect and export globals, but that detection has to follow calls several layers deep and fails quietly when it misses one.

Four arguments to `future_lapply` are load-bearing:

- `future.seed = TRUE` gives each task its own L'Ecuyer stream. Results are reproducible across runs and independent of worker count. **They do not match a sequential single-seed run, and cannot** — see the validation section.
- `future.chunk.size = 1` dispatches one split at a time. The default splits the range into equal blocks, which is wrong here because cost grows with the training window and the worker holding the last block finishes long after the others.
- `future.packages` lists the libraries explicitly. Workers inherit exported globals but not `library()` calls.
- `RhpcBLASctl::blas_set_num_threads(1)` inside the worker avoids oversubscription when several R processes each spawn linear-algebra threads.

**Checkpointing.** Each split writes its result to `_ckpt*/`, and a rerun loads completed splits instead of recomputing them. Errors are trapped per split, so one failure reports itself and leaves the rest intact rather than cancelling the run. Only successful results are checkpointed, so a failed split retries on the next run. Delete the checkpoint directory between genuinely different configurations or you will silently resume from the previous one.

------------------------------------------------------------------------

## 3. Validation sequence

The order matters. Changing the loop structure and the random-number handling at the same time makes any difference unattributable, so they were changed and checked separately.

**Stage 1. Refactor, still sequential.** The loop body became a function called with `lapply`. `lapply` consumes the random stream in the same order a `for` loop does, so this version must reproduce the original *bit for bit*. Checked with `identical()` on the forecast arrays over a short range. This step that establishes the restructuring didn't change the computation.

**Stage 2. Swap in `future_lapply`.** Per-task seeding necessarily changes the draws, so correctness can no longer be checked by comparison with the sequential run. Instead:

- run twice with identical settings, confirm `identical()` output (reproducible)
- run again with a different worker count, confirm `identical()` output (independent of how many cores happen to be free)

**Stage 3. Full run against the published tables.** Because the numbers can't match exactly, the criterion is that the parallel-versus-published difference is no larger than the difference between two legitimate runs under different seeds. Both exercises pass: roughly 0.5% median per cell either way, with interval coverage never moving more than 0.02.

**Stage 4. Restore MCMC settings.** Chain length was chosen from a mixing diagnostic (`dgfc_ess_diagnostic.R`) rather than by guessing a multiple of the old settings. At 5000 retained after 5000 burn-in, the latent correlations give near-independent draws and the margin adjustment gives about 200 effective draws.

**Stage 5. Factor sweep.** With the loops fast enough to run several times, the factor dimension was swept across values with the competitors switched off, since they do not depend on `k`.

------------------------------------------------------------------------

## 4. Files

### Modified helpers

**`_helpers/var.R`** — `stationary_var1_covariance` now symmetrizes its result before returning. The matrix is symmetric by construction, but the Kronecker solve loses that to floating-point rounding, and KFAS then fails its LDL decomposition when the result is passed as the initial state covariance `P1`. This surfaces at larger `k`, where the Kronecker system is bigger and worse conditioned. **Nothing runs above roughly `k = 5` without this fix.**

### New scripts

| File | Purpose |
|------------------------------------|------------------------------------|
| `application_forecasting_crime_parallel.R` | Parallel crime forecasting exercise. Produces Table 1. |
| `crime_v2_sweep.R` | Crime factor sweep. Takes `k` as a command-line argument, DGFC only. |
| `application_forecasting_macro_parallel.R` | Parallel macro forecasting exercise. Produces Table 2. |
| `macro_v2_sweep.R` | Macro factor sweep. Takes `k` as a command-line argument, DGFC only. |
| `dgfc_ess_diagnostic.R` | Mixing diagnostic. One long fit, effective sample sizes and trace plots for the identified quantities. |
| `dglm_diagnostic.R` | Importance-weight diagnostic for the Poisson and negative binomial competitors. |

Sweep scripts are run once per value:

``` bash
for k in 3 5 7 10; do Rscript crime_v2_sweep.R $k; done
for k in 3 7 10 14; do Rscript macro_v2_sweep.R  $k; done
```

### Environment

`scoringutils` must be pinned to **1.2.2**. Version 2.x removed `interval_score`, which `_helpers/forecasting.R` calls, so a fresh install breaks the scoring stage.

``` r
remotes::install_version("scoringutils", version = "1.2.2",
                         repos = "https://cloud.r-project.org")
```

`future`, `future.apply` and `RhpcBLASctl` are new dependencies.

------------------------------------------------------------------------

## 5. Results manifest (found in \*\_results folders)

Which settings produced which file.

### Crime

| File | Script | k | DGFC draw/burn | BVAR draw/burn | Seed | Runtime |
|-------------|----------------|------|----------|----------|---------|-----------|
| `crime_fcasts_og_repro.csv` | `crime_v2` (original) | 3 | 1000/0 | 1000/1000 | 8675309 | \~90 min seq |
| `crime_fcasts_par_stage2_*.csv` | `application_forecasting_crime_parallel` | 3 | 1000/0 | 1000/1000 | 8675309, 42065 | 10 min each |
| `crime_fcasts_stage4_*.csv` | `application_forecasting_crime_parallel` | 7 | 5000/5000 | 5000/5000 | 8675309, 42065 | 68 min each |
| `crime_fcasts_k{3,5,7,10}.csv` | `crime_v2_sweep` | 3/5/7/10 | 5000/5000 | off | 8675309 | 59, 60, 63, 68 min |

### Macro

| File | Script | k | DGFC draw/burn | BVAR draw/burn | Seed | Runtime |
|--------------|-------------|---------|----------|----------|---------|----------|
| `macro_og_repo.csv` | `macro_v1` (original) | 3 | 2000/0 | 1000/1000 | 8675309 | \~2.9 hr seq |
| `macro_table_full_V2_short_*.csv` | `application_forecasting_macro_parallel.R` | 3 | 2000/0 | 1000/1000 | 8675309, 42065 | 13.5 min each |
| `macro_table_stage4.csv` | `application_forecasting_macro_parallel.R` | 10 | 5000/5000 | 5000/5000 | 8675309 | 74 min |
| `macro_fcasts_k{3,7,10,14}.csv` | `macro_v2_sweep` | 3/7/10/14 | 5000/5000 | off | 8675309 | 57, 59, 61, 70 min |

The DLM and the two DGLM competitors are not MCMC: the DLM fits its variance parameters by maximum likelihood and draws from a single normal, and the DGLMs fit by maximum likelihood and draw by importance sampling. Their draw counts stayed at the original values throughout (2000 for the DLM, 1000 for each DGLM).

All runs on a 14-core Apple Silicon machine, 13 workers.

------------------------------------------------------------------------

## 6. Notes and gotchas

- Checkpoint directories (`_ckpt*/`), workspace images (`*.RData`) and logs are gitignored. Workspace images run to several hundred MB.
- `save.image` is called after the loop and before scoring in every script. A failure in the scoring stage would otherwise discard the whole run.
- The crime scoring window starts at `Tstart + H`, so shortening the loop below that for a test run makes `period` count backwards and silently scores the wrong rows. Macro has no equivalent trap: shortening `Fend` is consistent throughout.
- Per-split timings are printed from inside the workers, so they sum to far more than the wall clock and `message` output arrives in batches rather than streaming.

Ethan Carlson
