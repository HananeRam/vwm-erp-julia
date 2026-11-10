# Analysis settings, aligned with the MATLAB/EEGLAB/ERPLAB scripts (#0–#8a).
# Edit here, not in the pipeline code.

# --- Data ------------------------------------------------------------------
# Set the DS006979 environment variable, or put the dataset in data/ds006979.
const BIDS_ROOT = get(ENV, "DS006979", joinpath(@__DIR__, "..", "data", "ds006979"))
const TASK = "visualworkingmemorytask"

# --- Participants ----------------------------------------------------------
# All 53 participants (SUB list in Scripts #5, #7a, #8a).
const ALL_IDS = [1, 2, 3, 4, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                 16, 17, 18, 19, 21, 22, 23, 24, 25, 28, 30,
                 32, 33, 34, 35, 36, 39, 40, 41, 43, 46, 47, 48,
                 49, 50, 51, 52, 53, 54, 58, 59, 61, 62, 64,
                 69, 70, 72, 73, 74, 76]
# Final stimulus-locked sample (N = 39; same list as the alpha and NSW scripts).
const FINAL_39_IDS = [1, 2, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19,
                      21, 22, 23, 24, 28, 30, 32, 34, 35, 40, 47, 48, 49,
                      50, 51, 54, 58, 59, 61, 62, 64, 69, 70, 73, 76]
const SUBJECT_IDS = FINAL_39_IDS            # or ALL_IDS for all 53
const SUBJECTS = ["sub-" * lpad(string(i), 2, '0') for i in SUBJECT_IDS]

# Script #0 swapped response codes 120/121 for participants 21 and 22. The BIDS
# files were exported after that recoding, so nothing to swap here.
const SWAP_RESPONSE_IDS = Int[]

# --- Filters (Script #1 and #7a) --------------------------------------------
const HIGHPASS = (cutoff=0.1, order=2)      # Hz; ERPLAB order (12 dB/oct)
const LOWPASS  = (cutoff=20.0, order=8)     # Hz; ERPLAB order (48 dB/oct)
# As in Script #7a, the 20 Hz low-pass is applied to the epochs (equivalent to
# filtering the averaged ERPs), not to the continuous data.

# --- Sampling rate (Script #1) ----------------------------------------------
# As in Script #1 (pop_resample(EEG, 256)), every recording is downsampled to
# 256 Hz on loading, with an anti-aliasing filter. This also gives every
# participant the same time points (in dataset versions 1.0.x, the file of sub-49
# was stored at 2048 Hz). The measurement windows below (238.3 ms = sample 61,
# 996.1 ms = sample 255 at 256 Hz) are 256 Hz sample times.
const TARGET_FS = 256.0                     # Hz

# --- Epochs (Script #5) ------------------------------------------------------
const EPOCH    = (tmin=-0.2, tmax=1.0)      # seconds, locked to the memory array (codes 51–58)
const BASELINE = (-0.2, 0.0)

# --- Cluster of interest (Script #7a) ----------------------------------------
const CLUSTER = [("P3", "P4"), ("P5", "P6"), ("P7", "P8"), ("P9", "P10"),
                 ("PO7", "PO8"), ("PO3", "PO4"), ("O1", "O2")]

# Only these channels are read from disk (saves memory).
const LOAD_CHANNELS = unique(vcat([[l, r] for (l, r) in CLUSTER]...,
                                  ["HEOG-left", "HEOG-right", "VEOG-lower", "Fp2",
                                   "Mastoid-left", "Mastoid-right"]))

# --- Measurement windows (from the paper) --------------------------------------
# A posterior cluster spanned ~238–996 ms and was split into N2pc and CDA windows.
const N2PC_WINDOW = (0.2383, 0.300)         # seconds
const CDA_WINDOW  = (0.300, 0.9961)         # seconds, sustained CDA
# Early CDA (exploratory in the paper): transition from item selection to stable
# maintenance; linked to encoding and filtering efficiency.
const CDA_EARLY_WINDOW = (0.300, 0.450)     # seconds

# All six pairwise comparisons (A − B), in the order of the original R output
# (glht Tukey contrasts, Bonferroni-adjusted).
const ERP_PAIRS = [("3-RSP", "2-UR"), ("3-RWP", "2-UR"), ("3-UR", "2-UR"),
                   ("3-RWP", "3-RSP"), ("3-UR", "3-RSP"), ("3-UR", "3-RWP")]

# Reference values for the comparison table: the MATLAB/R analysis
# (correct trials only), as reported in the Stage 2 manuscript
# (2 decimals). Pairs and statistics the manuscript does not report are left out
# or NaN and print as "not reported" / "—".
# Mean, SD (µV) per condition (sustained CDA means are not reported).
const PAPER_MEANS = Dict(
    "n2pc"      => Dict("2-UR" => (-0.33, 0.67), "3-RSP" => (-0.66, 0.71),
                        "3-RWP" => (-0.57, 0.75), "3-UR" => (-0.33, 0.74)),
    "cda_early" => Dict("2-UR" => (-0.43, 0.78), "3-RSP" => (-0.43, 0.66),
                        "3-RWP" => (-0.73, 0.81), "3-UR" => (-0.52, 0.78)))
# Repeated-measures ANOVAs; p as reported.
const PAPER_ANOVA = Dict("n2pc"      => (F=7.32, df1=3, df2=114, p="< .001", eta_p2=0.16),
                         "cda_early" => (F=5.27, df1=3, df2=114, p="= .002", eta_p2=0.122),
                         "cda"       => (F=0.86, df1=3, df2=114, p="= .46",  eta_p2=0.02))
# Pairwise (A, B) => (mean difference, CI low, CI high, Bonferroni p);
# NaN = not reported (only p given in the manuscript).
const PAPER_PAIRWISE = Dict(
    "n2pc" => Dict(("3-RSP", "2-UR")  => (-0.33, -0.50, -0.15, 0.001),
                   ("3-RWP", "2-UR")  => (-0.24, -0.41, -0.07, 0.037),
                   ("3-UR", "3-RSP")  => ( 0.33,  0.15,  0.50, 0.001),
                   ("3-UR", "3-RWP")  => ( 0.24,  0.07,  0.41, 0.038)),
    "cda_early" => Dict(("3-RSP", "2-UR")  => (NaN, NaN, NaN, 1.0),
                        ("3-RWP", "2-UR")  => (-0.29, -0.46, -0.12, 0.004),
                        ("3-UR", "2-UR")   => (NaN, NaN, NaN, 1.0),
                        ("3-RWP", "3-RSP") => (-0.29, -0.46, -0.13, 0.003),
                        ("3-UR", "3-RSP")  => (NaN, NaN, NaN, 1.0),
                        ("3-UR", "3-RWP")  => (NaN, NaN, NaN, 0.102)))

# --- Artifact detection (Script #6a) -----------------------------------------
# Values from the AR_*_parameters_*.xlsx spreadsheets. Channel numbers there map
# to: 70 = bipolar HEOG, 71 = bipolar VEOG, [21:27 58:64] = the 14 cluster channels.
# Per-participant thresholds (blink, saccade_N2pc, saccade_CDA, ...) are in
# artifact_overrides.csv (columns: id, name, threshold), generated from the
# same spreadsheets; values below are the defaults.
const ARTIFACTS = [
    (name="blink",          detector=:step,    channel="VEOG",   period=(-0.2, 0.2),
     threshold=50.0,  window_ms=200, step_ms=10),
    (name="saccade_N2pc",   detector=:step,    channel="HEOG",   period=(-0.05, 0.3),
     threshold=16.0,  window_ms=200, step_ms=10),
    (name="saccade_CDA",    detector=:step,    channel="HEOG",   period=(0.3, 1.0),
     threshold=16.0,  window_ms=200, step_ms=10),
    (name="saccade_epoch",  detector=:step,    channel="HEOG",   period=(-0.2, 0.95),
     threshold=64.0,  window_ms=200, step_ms=10),
    (name="crap_extreme",   detector=:extreme, channel=:cluster, period=(-0.2, 1.0),
     threshold=200.0, window_ms=0,   step_ms=0),
    (name="crap_mwpp",      detector=:mwpp,    channel=:cluster, period=(-0.2, 1.0),
     threshold=150.0, window_ms=200, step_ms=100),
]

# --- Blink correction (Scripts #2–#4) -------------------------------------------
#   :ica_weights  the original ICA solution: the EEGLAB weights in
#                 <id>_allERPs_optimized_weights.set (dataset derivatives, searched
#                 under ICA_WEIGHTS_DIR) and the blink components removed in
#                 ICs_to_Remove.csv. Reproduces Script #4.
#   :fastica      a new, automated ICA in Julia (removes as many components as
#                 the original, those most correlated with VEOG). Approximate.
#   :reject       no ICA: blinks after 200 ms rejected (30 µV, Diaz et al., 2021).
const BLINK_CORRECTION = :ica_weights
const ICA_WEIGHTS_DIR = get(ENV, "ICA_WEIGHTS",
                            joinpath(BIDS_ROOT, "derivatives", "ReplicationPeterson2015"))
const ICS_TO_REMOVE = joinpath(@__DIR__, "ICs_to_Remove.csv")
# Components actually removed by Script #4 where they differ from ICs_to_Remove.xlsx,
# identified by reconstructing the MATLAB ..._ica_corr_reref.set files from the
# uncorrected data and the weights (see vwm-erp-julia-probe). Julia removes these to
# reproduce the published analysis. (Participant 46, also affected, is not in this sample.)
const ICS_APPLIED = Dict(49 => [2], 69 => [1, 7])
const USE_ICA = BLINK_CORRECTION !== :reject
const ICA_N_REMOVE = Dict(3 => 2, 4 => 2, 36 => 2, 39 => 2, 40 => 2, 32 => 0,
                          20 => 0, 26 => 0, 27 => 0, 29 => 0, 31 => 0, 37 => 0, 38 => 0,
                          42 => 0, 44 => 0, 45 => 0, 55 => 0, 56 => 0, 57 => 0, 60 => 0,
                          63 => 0, 65 => 0, 66 => 0, 67 => 0, 68 => 0, 71 => 0, 75 => 0)
ica_n_remove(id) = get(ICA_N_REMOVE, id, 1)       # for :fastica
const ICA_COMPONENTS = 30                          # for :fastica

# Late blinks: used instead of ICA when USE_ICA = false. Same step function as
# the other checks (200 ms window, 10 ms steps) on VEOG after 200 ms, with the
# 30 µV threshold of Diaz, Vogel & Awh (2021). The −200 to 200 ms blink check
# above (blinks during stimulus presentation) is applied in both cases.
const USE_LATE_BLINK_REJECTION = !USE_ICA
const LATE_BLINK = (name="blink_late", detector=:step, channel="VEOG", period=(0.2, 1.0),
                    threshold=30.0, window_ms=200, step_ms=10)
USE_LATE_BLINK_REJECTION && push!(ARTIFACTS, LATE_BLINK)
const ARTIFACT_OVERRIDES = joinpath(@__DIR__, "artifact_overrides.csv")

# BDF_stimulus_correct.txt bins require a correct response (120) after the memory array.
const CORRECT_ONLY = true

const CONDITIONS = ["2-UR", "3-UR", "3-RSP", "3-RWP"]   # order used in the figures
# Same colors as the MATLAB figures (MATLAB default line colors).
const COLORS = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E"]

# --- Exploratory analyses (MATLAB Alpha and NSW scripts) -------------------------
# Both use the non-lateralized average of the 14 cluster channels, from the
# artifact-free epochs before the 20 Hz low-pass.
const ALPHA_BAND     = (8.0, 12.0)          # Hz
const ALPHA_PAD_S    = 0.5                  # mirror padding before filtering (s)
# false = as in the MATLAB script: Hilbert transform after the padding is
# removed, which creates edge artifacts in the first and last ~50 ms (inside
# the baseline and the end of the 300–996.1 ms window). true = Hilbert
# transform on the padded signal, which avoids them (sensitivity analysis).
const ALPHA_HILBERT_PADDED = false
# Band-pass design:
#   :matlab      the exact filter of MATLAB's bandpass(x, [8 12], 256), read from
#                scripts/matlab_bandpass_coeffs.csv (exported in MATLAB, see README)
#   :elliptic    an approximation of bandpass() defaults
#   :butterworth a gentle filter
const ALPHA_MATLAB_COEFFS = joinpath(@__DIR__, "matlab_bandpass_coeffs.csv")
const ALPHA_FILTER = isfile(ALPHA_MATLAB_COEFFS) ? :matlab : :elliptic
# Candidate designs computed as well (alpha_<design>.csv), to find which one
# MATLAB's bandpass() used, with scripts/05_identify_alpha_filter.jl.
const ALPHA_CANDIDATES = [:elliptic, :cheby1, :cheby2, :butter15, :butterworth]
# MATLAB alpha results (Step1_Compute_Alpha_Power output), for the comparison.
const ALPHA_MATLAB_MAT = get(ENV, "ALPHA_MAT", joinpath(@__DIR__, "alpha_power.mat"))
const ALPHA_WINDOW   = (0.300, 0.9961)      # seconds, matches the CDA window
const PAIRWISE = [("2-UR", "3-UR"), ("2-UR", "3-RSP"), ("2-UR", "3-RWP"),
                  ("3-UR", "3-RSP"), ("3-UR", "3-RWP")]

# Alpha results (correct trials): BF01 from BayesFactor::ttestBF
# (paired) in R, on % change averaged over 300–996.1 ms.
const PAPER_ALPHA_BF01 = Dict(("2-UR", "3-UR")  => 5.49, ("2-UR", "3-RSP") => 5.79,
                              ("2-UR", "3-RWP") => 1.00, ("3-UR", "3-RSP") => 5.27,
                              ("3-UR", "3-RWP") => 2.82)

const NSW_PERMUTATIONS = 10_000
const NSW_ALPHA        = 0.05
const NSW_FIGURE_LOWPASS = 30.0             # Hz, figure only
# Windows from the cluster tests (Step 4, correct trials) and the
# comparisons tested in each.
const NSW_WINDOWS = Dict("Grouping_Cluster"    => (0.3281, 0.6172),
                         "SetSize_Cluster"     => (0.4023, 0.6055),
                         "Late_Cluster"        => (0.8867, 0.9609),
                         "Diaz_SetSize"        => (0.330, 0.880),
                         "SetSize_Cluster_RWP" => (0.3242, 0.5742))
const NSW_TARGETS = [("Grouping_Cluster", "3-UR", "3-RSP"),
                     ("Grouping_Cluster", "3-UR", "3-RWP"),
                     ("SetSize_Cluster", "2-UR", "3-UR"),
                     ("SetSize_Cluster_RWP", "2-UR", "3-RWP"),
                     ("Late_Cluster", "2-UR", "3-UR"),
                     ("Diaz_SetSize", "2-UR", "3-RSP")]
# NSW results (correct trials). Significant clusters (ms) per comparison
# (none = no cluster), and BF10 from BayesFactor::ttestBF (paired) for each window
# test, in NSW_TARGETS order (BF01 = 5.22 and 2.99 in the manuscript → BF10 = 1/BF01).
const PAPER_NSW_CLUSTERS = Dict(("2-UR", "3-UR")  => [(402.3, 605.5), (886.7, 960.9)],
                                ("2-UR", "3-RSP") => Tuple{Float64,Float64}[],
                                ("2-UR", "3-RWP") => [(324.2, 574.2)],
                                ("3-UR", "3-RSP") => [(328.1, 617.2)],
                                ("3-UR", "3-RWP") => Tuple{Float64,Float64}[])
const PAPER_NSW_BF10 = Dict(("Grouping_Cluster", "3-UR", "3-RSP")    => 23892.0,
                            ("Grouping_Cluster", "3-UR", "3-RWP")    => 1 / 5.22,
                            ("SetSize_Cluster", "2-UR", "3-UR")      => 102.66,
                            ("SetSize_Cluster_RWP", "2-UR", "3-RWP") => 27.91,
                            ("Late_Cluster", "2-UR", "3-UR")         => 5.09,
                            ("Diaz_SetSize", "2-UR", "3-RSP")        => 1 / 2.99)
const RESULTS = joinpath(@__DIR__, "..", "results")
