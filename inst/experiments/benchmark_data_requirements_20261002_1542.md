# Data wanted for benchmarking batchmix

## Purpose

We are looking for datasets on which to benchmark `batchmix`, a Bayesian mixture model that classifies items from multivariate measurements while correcting for batch effects. This note says what the model needs, which properties make a dataset informative, and what we have already tried. The data do not have to be ELISA. Any item-level, multivariate measurement with a recorded batch and a latent discrete class will do.

## What the model does

Each item has a vector of measurements, a batch label, and an unobserved class (K classes; typically negative/positive, but any small number).

- **Classes.** Within a class the measurements are multivariate normal (or Student-t) on a scale the user chooses, usually log. Some items can have known labels, which stay fixed.
- **Batch effects.** Each batch has an additive shift on the class means and a multiplicative inflation of the variance. An optional batch-by-class interaction lets the shift differ by class. Shifts and scales are partially pooled across batches.
- **Class proportions per batch.** One of: shared across batches (global), partially pooled (exchangeable batches), or smooth over a one-dimensional coordinate such as collection date (Gaussian process, default a random walk; non-stationary by design).
- **Optional.** Binary (probit) and left- or right-censored columns are supported (`MVN_MIXED`); missing values are treated as ignorable.
- **Outputs.** Class allocations, batch-corrected data, per-batch class proportions, and predictions for a new batch.

Limits to keep in mind: one batch variable per fit (so a technical factor and an epidemiological factor cannot yet be modelled separately); Gaussian class shapes, so skewed or heavy-tailed classes degrade the fit; no covariates on the class proportions.

## What the data must contain

| Requirement | Why |
|---|---|
| **Item-level measurements, raw**, not aggregated counts and not already batch-normalised | Batch correction is the point. Plate-adjusted values hide the batch signal we want to estimate. |
| **At least two correlated measurements per item** | A single feature gave no separate positive mode in our SCAPES check, so the class split rested on the prior. |
| **A batch identifier on every item** (plate, run, instrument, lab, site, reagent lot, assay date) | The model's batch terms are defined by it. |
| **Known class truth for a subset of items, from an independent reference method** (PCR, culture, neutralisation, clinical diagnosis, a mixed-at-known-ratio panel) | Without this we can only compare models with each other. |
| **Both classes present in the same batches** | If negatives and positives sit on different plates, class and batch effects cannot be separated. The Kenyan InBios/KWTRP comparison has this flaw. |
| **At least 10 batches** (preferably 15 or more), **at least about 30 items per batch** | Pooled batch hyperparameters and the weight hierarchy need enough batches to estimate a spread. These thresholds are our judgement from the experiments so far. We have not tested the lower limit. |
| **Open licence and no access gate** (for example CC BY, direct download) | Several candidate sets needed a guestbook or access request. |

## Strongly preferred

1. **Replicate or bridging samples across batches**, with the batch ID on each replicate: the same sample or the same control material on every plate or run. This is the cleanest check on batch correction, because replicates of one sample should agree after correction. Stockholm's "Patient 4" (24 measurements of one person) is useful in this way but has no plate IDs, which limits it.
2. **Per-batch control wells or controls** (negative, positive, blank) with their raw readings.
3. **A time or space coordinate for every batch**, with at least 10 distinct values, irregularly spaced, and some real abrupt changes (an outbreak, a policy or reagent change) as well as gradual drift. This tests the Gaussian process weight prior.
4. **Per-batch class proportions that are known or independently verified**, for example because every item has a reference-method result, or because we can build semi-synthetic batches from gold-standard items at chosen proportions.
5. **Class separation that is neither trivial nor hopeless.** As a rough guide, an AUC between about 0.85 and 0.98 for the best single feature. With near-perfect separation (Stockholm, resampled positive controls) every weight prior gave the same answer, so the data cannot discriminate between our models. With one weak feature we cannot identify the classes.
6. **Batch effects that matter**: a between-batch shift comparable to the within-class standard deviation, visible in the controls.
7. **A batch that can be held out**, ideally the most recent one, for forecasting and new-batch prediction.
8. **Both the technical and the epidemiological grouping recorded**, if they differ (plate and week, for instance). We cannot fit both at once yet, but we want the choice to be ours.

## Modalities that fit

Not verified as specific datasets; these are the kinds of data that match the design.

- **Multiplex or bead-based serology** (median fluorescence intensity for several antigens per sample, with plate and bead-lot IDs). This is probably the best fit: several correlated features, plates as batches, and often a gold-standard panel.
- **Flow or mass cytometry**: marker intensities per cell or per cell cluster, batch = staining or acquisition run, bridge samples included, class = a gated or manually labelled population.
- **qPCR or multiplex PCR** with Ct values: right-censored at the cycle limit, which uses the censoring support. Batch = run or instrument.
- **Clinical laboratory analytes measured at several sites or analysers**, class = diagnosis, batch = lab or analyser.
- **Mixed binary and continuous features** (for example rapid-test readouts plus quantitative values).
- **Single-cell or omics summaries** (principal components) with batch and a labelled cell type, if the classes are approximately Gaussian after transformation.

A single-feature assay is not enough on its own, however good the labels.

## What we have already tried, and what it told us

| Dataset | Result | Why it falls short |
|---|---|---|
| Stockholm blood donors (raw spike and RBD OD, 21 weeks, controls) | Held-out control sensitivity 0.986 and specificity about 0.98, and all weight priors agreed. | Signal too strong to separate priors; weeks nested in four assay runs; no plate IDs; truth only from controls. |
| Stockholm, semi-synthetic from the controls | All priors tied (absolute error 0.002 to 0.006 in weekly prevalence). | Same reason. |
| SCAPES Lassa IgG (33 plates) | Not fitted. | One feature, no gold standard, no control wells, plate and time crossed but the model cannot split them. |
| Kenya InBios/KWTRP comparison | Not fitted. Has gold-standard negatives and positives and two assays. | Classes sit on different plates; no time structure. |
| Malawi (counts only) | Not fitted. | No item-level values; no truth. |

## How we will evaluate on a new dataset

- Held-out labelled items: sensitivity, specificity, log score, calibration.
- Batch correction: spread of replicates or controls before and after correction.
- Per-batch class proportions against the reference values: absolute error and interval coverage.
- Posterior predictive checks: per-batch means, standard deviations and upper tails.
- Forecasting a held-out batch (time-smoothed weights only).
- Convergence (Rhat, effective sample size) and sensitivity to prior choices. We currently see non-convergence when the interaction term is on, which new data should help us diagnose.

## Format we would like

One row per item, as CSV, with a short data dictionary:

- `item_id`; `batch_id`; optional `batch_coordinate` (date, or latitude/longitude).
- The raw feature columns, with units and any detection limits; a flag for censored values.
- `class_label`, with `NA` where unknown, and `label_source` (which reference method).
- `replicate_group` for repeats of the same sample or control; `well_type` (sample, positive control, negative control, blank).
- Optional covariates (age, sex, site).

Please also send the study protocol or paper, the licence, and any known issues (duplicated IDs, reruns, reagent changes).

## Quick checklist for a candidate dataset

- [ ] Item-level raw values, at least two features
- [ ] Batch ID on every item, at least 10 batches
- [ ] Reference-method truth for a subset, with both classes in most batches
- [ ] Replicates or controls with batch IDs
- [ ] Time or space coordinate per batch (if testing smoothing)
- [ ] Open licence, direct download
- [ ] Not near-perfectly separable, not a single weak feature
