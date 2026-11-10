# Lateralized ERPs in visual working memory: a Julia reanalysis

A small, reproducible Julia pipeline that reanalyzes my open EEG dataset
[ds006979 on OpenNeuro](https://doi.org/10.18112/openneuro.ds006979.v1.1.4) (version 1.1.4)
(Ramzaoui & Beck, *Scientific Reports*, Registered Report; Stage 1 accepted, Stage 2 under review).

The study used a bilateral change-detection task to test how perceptual grouping by
color repetition affects stages of visual working memory. This repository computes two
lateralized ERP markers from the raw EEG:

- **N2pc**: individuation during encoding
- **CDA** (contralateral delay activity): maintenance and memory load

These are the **stimulus-locked** analyses (time-locked to the memory array). The
**probe-locked ERPs** of the original study (**N2pc, N2, FN400**) are reproduced in a separate
repository, [vwm-erp-julia-probe](https://github.com/HananeRam/vwm-erp-julia-probe).

Two exploratory analyses are also reproduced, on the non-lateralized average of the
14 posterior channels:

- **Alpha-power suppression**: 8–12 Hz Hilbert power, % change from the −200–0 ms baseline,
  averaged over 300–996.1 ms; Bayesian paired t-tests (BF01, Cauchy prior r = √2/2, computed in
  Julia exactly as `BayesFactor::ttestBF` in R), plus paired t-tests with Bonferroni correction
- **Negative slow wave (NSW)**: cluster-based permutation tests (10,000 sign flips) on
  paired difference waves, plus mean amplitudes in the cluster windows with Bayesian paired
  t-tests (BF10/BF01)

Measurement windows follow the paper: **N2pc 238.3–300 ms**, **sustained CDA 300–996.1 ms**,
and **early CDA 300–450 ms** (exploratory in the paper; it isolates the transition from item
selection to stable maintenance and is linked to encoding and filtering efficiency)
(final sample N = 39).

> This is an independent demonstration of the analysis in Julia. The confirmatory
> analysis is the one specified in the Stage 1 protocol
> ([OSF 8ZS96](https://doi.org/10.17605/OSF.IO/8ZS96)), run in MATLAB (EEGLAB/ERPLAB).

## Pipeline

The steps follow the original MATLAB/EEGLAB/ERPLAB scripts used for the registered analysis.

| MATLAB script | Step | Julia |
| --- | --- | --- |
| #0 | BDF import, channel selection, CleanLine (60 Hz), practice-trial removal, response recoding (participants 21, 22) | Already applied: the BIDS files are the Script #0 output, so preprocessing starts at Script #1 |
| #1 | Downsampling to 256 Hz; mastoid reference; bipolar HEOG (left − right) and VEOG (lower − Fp2); DC removal; 0.1 Hz high-pass (Butterworth, 12 dB/oct) | `resample`, `rereference!`, `add_bipolars!`, `remove_dc!`, `butter_filter!` |
| #2–#4 | ICA correction of blinks; interpolation of bad channels | **Original EEGLAB ICA weights** from the OpenNeuro derivatives and the removed components in `scripts/ICs_to_Remove.csv` (`read_ica_weights`, `apply_ica_weights!`). Interpolation not needed: no channel of interest was interpolated in the original analysis |
| #5 | Trials binned by condition from event codes; epochs −200 to 1000 ms; baseline −200 to 0 ms | `label_trials`, `epoch` |
| #6a | Blink and saccade step functions (VEOG, HEOG), extreme values, moving-window peak-to-peak | `art_step`, `art_extreme`, `art_mwpp` |
| #7a | Contra − ipsi waves averaged over P3/4, P5/6, P7/8, P9/10, PO7/8, PO3/4, O1/2; left-cue and right-cue trials averaged separately and combined 50/50, as in the ERPLAB bin operations (`nbin = 0.5*bin_left + 0.5*bin_right`); 20 Hz low-pass (48 dB/oct) on the averages | `cluster_contra_minus_ipsi`, `lowpass_epochs!` |
| #11a | Mean amplitudes (`meanbl`, baseline −200–0 ms) on the averages **without** the low-pass; window edges at the closest samples | `erplab_window`, `results/participant_means.csv` |
| #8a | HEOG difference waves to check residual eye movements | `results/heog_check.png` |
| Paper | One-way repeated-measures ANOVAs and pairwise comparisons on participant means (N2pc, early CDA 300–450 ms, sustained CDA); side-by-side comparison with the MATLAB values | `rm_anova`, `pairwise_pooled`, `results/paper_comparison.txt` |
| New | Single-trial mixed models; model-based effects | MixedModels.jl, [Effects.jl](https://github.com/beacon-biosignals/Effects.jl) |
| Alpha Steps 1–3 | ROI alpha power (8–12 Hz, mirror padding, Hilbert), % change, paired t-tests, exports for Bayesian tests | `alpha_power`, `src/ExploratoryStats.jl`, `scripts/04_…` |
| NSW Steps 1–4 | ROI slow wave; cluster permutation tests; window extraction and paired tests | `cluster_permutation`, `scripts/04_…` |

**Data source.** EEG was recorded with a BioSemi ActiveTwo system (64 scalp channels, 3 EOG,
2 mastoids; 512 Hz) and saved as BioSemi `.bdf` files. Script #0 imported the `.bdf` files in
EEGLAB (channel selection, CleanLine, practice-trial removal), and the cleaned recordings were
converted to BIDS with FieldTrip. The BIDS files therefore correspond to the Script #0
output (`N.set`), including the recoded responses for participants 21 and 22, and
preprocessing starts at Script #1, as in the original analysis. The BIDS release provides them as BrainVision files for every
participant (and also as EDF for sub-01). This repository reads the BrainVision files with a small
memory-mapped reader; `scripts/03_check_formats.jl` reads sub-01's EDF copy with
[EDF.jl](https://github.com/beacon-biosignals/EDF.jl) and confirms both formats give the same
signals (r = 1.000 on every channel; differences ≤ 0.55 µV, from EDF's 16-bit storage).

All settings (participants, filters, cluster, windows, artifact thresholds) are in `scripts/config.jl`.

## Get the data

Every participant has BrainVision files (`.vhdr`, `.vmrk`, `.eeg`, about 600–720 MB each);
sub-01 also has an EDF copy. `scripts/download_data.sh` downloads the 39 participants of the final
sample into `data/ds006979` (about 25 GB; files already downloaded are skipped, so it can be rerun
after an interruption):

```bash
bash scripts/download_data.sh
```

The ICA weights of the original analysis are in the dataset derivatives
(`derivatives/ReplicationPeterson2015/sub-XX/<id>_allERPs_optimized_weights.set`, where `<id>` is
the participant number without leading zero, e.g. `sub-01/1_allERPs_optimized_weights.set`).
Download them for the same 39 participants (about 105 MB each):

```bash
for s in 01 02 07 08 09 10 11 12 13 14 16 17 18 19 21 22 23 24 28 30 32 34 35 40 \
         47 48 49 50 51 54 58 59 61 62 64 69 70 73 76; do
  n=$((10#$s)); d=data/ds006979/derivatives/ReplicationPeterson2015/sub-$s; mkdir -p $d
  curl -sL -o "$d/${n}_allERPs_optimized_weights.set" \
    https://s3.amazonaws.com/openneuro.org/ds006979/derivatives/ReplicationPeterson2015/sub-$s/${n}_allERPs_optimized_weights.set
done
```

The components removed for each participant (`ICs_to_Remove.xlsx`, also in the derivatives) are
included here as `scripts/ICs_to_Remove.csv`.

Or download from the [OpenNeuro page](https://openneuro.org/datasets/ds006979).
To use a dataset stored elsewhere, set `DS006979=/path/to/ds006979` (and, if the ICA weights are
elsewhere too, `ICA_WEIGHTS=/path/to/ReplicationPeterson2015`).

## Run

Requires [Julia](https://julialang.org/downloads/) 1.10 or later.

```bash
julia --project=. scripts/00_setup.jl                 # install packages (once)
julia --project=. scripts/01_process_subjects.jl      # preprocessing → results/*.csv
julia --project=. scripts/02_statistics_and_figures.jl # models and figures
julia --project=. scripts/03_check_formats.jl          # EDF vs BrainVision check (sub-01)
julia --project=. scripts/04_exploratory_alpha_nsw.jl  # alpha power and negative slow wave
julia --project=. scripts/05_identify_alpha_filter.jl  # check the alpha filter against MATLAB (full comparison needs your alpha_power.mat)
```

## Outputs (`results/`)

| File | Contents |
| --- | --- |
| `qc.csv` | Binned (correct) trials, rejection rate per detector and overall, per participant |
| `trials.csv` | Every trial: labels, accuracy, one flag per artifact detector, kept or not, N2pc and CDA |
| `participant_means.csv`, `condition_summary.csv` | Participant-level mean amplitudes and condition summaries |
| `heog.csv`, `heog_check.png` | HEOG contra − ipsi waves (eye-movement check) |
| `rejection_rates.png` | Percentage of trials rejected per participant |
| `waveforms.csv` | Participant × condition average lateralized waveforms |
| `erp_anova.csv`, `erp_pairwise.csv` | Repeated-measures ANOVAs and pairwise comparisons (N2pc, early and sustained CDA), same method as the R analysis |
| `erp_pairwise_paired_t.csv` | The same comparisons as paired t-tests (Bonferroni and Holm), for reference |
| `paper_comparison.txt` | Julia results next to the MATLAB values |
| `models.txt` | Mixed-model summaries for CDA and N2pc |
| `cda_effects.csv`, `cda_effects.png` | Model-based CDA by condition (Effects.jl) |
| `grand_average.png`, `grand_average_early.png` | Grand-average contra − ipsi waveforms by condition, in the style of the paper's Figures 4 (N2pc and sustained CDA) and 3 (N2pc and early CDA, −200 to 500 ms) |
| `exploratory/` | Alpha and NSW tables (t-tests and Bayes factors), `alpha_paper_comparison.txt`, `nsw_paper_comparison.txt`, cluster results, figures, and paired-data CSVs for Bayesian t-tests in R |

## Figures

Grand-average contra − ipsi waves in the posterior cluster, with the N2pc and sustained CDA
windows (`results/grand_average.png`; the N2pc and early CDA version is `grand_average_early.png`):

![Grand-average contra-minus-ipsi waves](results/grand_average.png)

Time-resolved alpha power (8–12 Hz, % change from baseline) and mean suppression in
300–996 ms (`results/exploratory/alpha_figure.png`):

![Alpha power](results/exploratory/alpha_figure.png)

Negative slow wave: grand averages with the MATLAB cluster windows, and mean amplitudes in
the primary maintenance and late set-size clusters (`results/exploratory/nsw_figure.png`):

![Negative slow wave](results/exploratory/nsw_figure.png)

## Validation against the MATLAB analysis (N = 39)

Applied to the BIDS data (version 1.1.4) with the original ICA weights, the Julia pipeline
reproduces the MATLAB/ERPLAB + R analysis (correct trials only). The reference is the revised
MATLAB analysis, in which the blink components removed for every participant are those listed in
`ICs_to_Remove.xlsx`. Julia keeps the same 31,954 trials as MATLAB, every condition
mean is identical (to 2 decimals), and every ANOVA, pairwise comparison, cluster after 300 ms
and Bayes factor gives the MATLAB/R value. `results/paper_comparison.txt` and
`results/exploratory/*_paper_comparison.txt` print both side by side.

| Analysis | MATLAB/R | Julia | Reproduced |
| --- | --- | --- | --- |
| Trials kept | 31,954 | 31,954 | Yes |
| N2pc, ANOVA | F(3, 114) = 7.44, p < .001, ηp² = .16 | F(3, 114) = 7.44, p < .001, ηp² = .16 | Yes |
| N2pc, 3-RSP vs 2-UR | −0.33 µV, p = .001 | −0.33 µV, p = .001 | Yes |
| N2pc, 3-RWP vs 2-UR | −0.24 µV, p = .038 | −0.24 µV, p = .038 | Yes |
| N2pc, 3-UR vs 3-RSP | 0.33 µV, p = .001 | 0.33 µV, p = .001 | Yes |
| N2pc, 3-UR vs 3-RWP | 0.24 µV, p = .037 | 0.24 µV, p = .037 | Yes |
| Early CDA, ANOVA | F(3, 114) = 5.37, p = .002, ηp² = .12 | F(3, 114) = 5.37, p = .002, ηp² = .12 | Yes |
| Early CDA, 3-RWP vs 2-UR and vs 3-RSP | −0.29 and −0.29 µV, p = .003 and .003 | −0.29 and −0.29 µV, p = .003 and .003 | Yes |
| Early CDA, 3-UR vs 3-RWP | 0.21 µV, p = .090 | 0.21 µV, p = .089 | Yes |
| Sustained CDA, ANOVA | F(3, 114) = 0.91, p = .44, ηp² = .02 | F(3, 114) = 0.90, p = .44, ηp² = .02 | Yes |
| NSW, grouping cluster (3-UR vs 3-RSP) | 328.1–617.2 ms, BF10 = 15,321 | 328.1–617.2 ms, BF10 = 15,304 | Yes |
| NSW, 3-UR vs 3-RWP | no cluster, BF01 = 5.42 | no cluster, BF01 = 5.42 | Yes |
| NSW, set size (2-UR vs 3-UR) | 402.3–605.5 ms, BF10 = 76.8 | 402.3–605.5 ms, BF10 = 76.8 | Yes |
| NSW, set size (2-UR vs 3-RWP) | 324.2–574.2 ms, BF10 = 31.3 | 324.2–574.2 ms, BF10 = 31.3 | Yes |
| NSW, late cluster (2-UR vs 3-UR) | 820.3–960.9 ms, BF10 = 7.21 | 820.3–960.9 ms, BF10 = 7.20 | Yes |
| NSW, 2-UR vs 3-RSP (a priori set-size window, 330–880 ms) | no cluster, BF01 = 3.11 | no cluster after 300 ms, BF01 = 3.11 | Yes |
| Alpha, BF01 (2-UR vs 3-UR, 3-RSP, 3-RWP; 3-UR vs 3-RSP, 3-RWP) | 5.49, 5.79, 0.95, 5.27, 2.74 | 5.49, 5.79, 0.95, 5.27, 2.74 | Yes |

Window tests use the MATLAB cluster windows, so both pipelines test the same intervals. The
last-digit differences (e.g. sustained CDA F = 0.905 vs 0.90, BF10 = 15,321 vs 15,304) come from
the MATLAB measurement and export files, which store values with 3 decimals, while Julia keeps
full precision; the trials, condition means and cluster windows are identical.

**Dataset version.** Use the latest version of the dataset, 1.1.4: it contains all the
corrections to the dataset and its files.

Julia also finds short early NSW clusters (about 130–300 ms) in the comparisons with 2-UR. The
original analysis found them too but did not report them, because they reflect sensory
differences between two- and three-item displays rather than working memory load. The table
above lists only the clusters after 300 ms.

## Notes

- **Probe-locked ERPs.** The probe-locked analyses of the study (N2pc, N2, FN400) are reproduced in
  [vwm-erp-julia-probe](https://github.com/HananeRam/vwm-erp-julia-probe) (N = 53).

## Data and citation

Data: Ramzaoui, H., & Beck, M. (2026). *Examining perceptual grouping on stages of
processing in visual working memory: An ERP study* [Dataset]. OpenNeuro, ds006979,
v1.1.4. https://doi.org/10.18112/openneuro.ds006979.v1.1.4 (CC0)

Code: MIT License.
