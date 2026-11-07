# Alpha band-pass filter: check against MATLAB.
# Run after 01:  julia --project=. scripts/05_identify_alpha_filter.jl
#
# 1. Always: checks that Julia's `matlab_filtfilt` with MATLAB's coefficients
#    (scripts/matlab_bandpass_coeffs.csv) reproduces MATLAB's filtfilt on a MATLAB
#    test signal (scripts/matlab_filtfilt_test.csv).
# 2. Only with your MATLAB alpha output (alpha_power.mat from Step 1, not part of
#    the repository): computes alpha power on the ROI trials saved by script 01
#    (results/roi_trials.jls) with MATLAB's filter and with candidate filters, and
#    compares each with the MATLAB output participant by participant.
#    Put alpha_power.mat in scripts/ or set ALPHA_MAT=/path/to/alpha_power.mat.

include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "..", "src", "VWMPipeline.jl"))
using .VWMPipeline
using CSV, DataFrames, Statistics, MAT, Printf, DSP, Serialization

# --- 1. filtfilt check on MATLAB's test signal ------------------------------------
test_file = joinpath(@__DIR__, "matlab_filtfilt_test.csv")
if isfile(ALPHA_MATLAB_COEFFS) && isfile(test_file)
    S_test = Matrix{Float64}(CSV.read(ALPHA_MATLAB_COEFFS, DataFrame; header=false))
    T = Matrix{Float64}(CSV.read(test_file, DataFrame; header=false))   # x, filtfilt(d, x)
    dmax = maximum(abs.(matlab_filtfilt(S_test, T[:, 1]) .- T[:, 2]))
    @printf("\nJulia's MATLAB-style filtfilt vs MATLAB on the test signal: max difference = %.2e (output range ±%.2f)\n",
            dmax, maximum(abs.(T[:, 2])))
end

# EEGLAB pop_resample check (512 → 256 Hz) on a MATLAB test signal
rin, rout = joinpath(@__DIR__, "eeglab_resample_in.csv"), joinpath(@__DIR__, "eeglab_resample_out.csv")
if isfile(rin) && isfile(rout)
    xin = vec(Matrix{Float64}(CSV.read(rin, DataFrame; header=false)))
    yref = vec(Matrix{Float64}(CSV.read(rout, DataFrame; header=false)))
    yj = eeglab_resample(xin, 2)
    if length(yj) == length(yref)
        @printf("Julia's eeglab_resample vs EEGLAB pop_resample: max difference = %.2e (signal range ±%.2f)\n",
                maximum(abs.(Float32.(yj) .- yref)), maximum(abs.(yref)))   # EEGLAB stores data in single precision
    else
        println("eeglab_resample: length $(length(yj)) vs EEGLAB $(length(yref))")
    end
end

# --- 2. Comparison with the MATLAB alpha output (optional) -------------------------
if !isfile(ALPHA_MATLAB_MAT)
    println("\nMATLAB alpha output not found ($ALPHA_MATLAB_MAT): skipping the comparison.")
    println("Put alpha_power.mat (saved by Step 1) in scripts/ or set ALPHA_MAT=/path/to/alpha_power.mat.")
    exit(0)
end
trials_file = joinpath(RESULTS, "roi_trials.jls")
isfile(trials_file) || error("$trials_file not found: rerun script 01")

# --- MATLAB output --------------------------------------------------------------
M = matread(ALPHA_MATLAB_MAT)
mtimes = vec(M["times"])                                   # ms
msubs = ["sub-" * lpad(strip(string(s)), 2, '0') for s in vec(M["SUB"])]
mgroups = [string(g) for g in vec(M["groupNames"])]
mpower = M["alphaPowerTimePerSubjectPerGroup"]
base = (mtimes .>= -200) .& (mtimes .<= 0)
win  = (mtimes .>= 300) .& (mtimes .<= 996.1)
pct(p) = 100 .* (p .- mean(p[base])) ./ mean(p[base])
mref = Dict((msubs[i], mgroups[g]) => vec(mpower[i, g]) for i in eachindex(msubs), g in eachindex(mgroups))

# --- ROI trials (Julia, with the original ICA) -----------------------------------
R = deserialize(trials_file)
fs = first(values(R)).fs
nyq = fs / 2
jt = 1000 .* first(values(R)).times
idx = [findfirst(t -> abs(t - m) < 0.01, jt) for m in mtimes]
any(isnothing, idx) && error("Time points differ between MATLAB and Julia (rerun script 01)")

# --- Candidate filters --------------------------------------------------------
# Passband 8–12 Hz; MATLAB bandpass() defaults: steepness 0.85 (stopbands below
# 6.8 and above 29.4 Hz at 256 Hz), 0.1 dB passband ripple, 60 dB stopband.
iir(proto, f1=8.0, f2=12.0) = convert(SecondOrderSections,
                                      digitalfilter(Bandpass(f1 / nyq, f2 / nyq), proto))
cands = Tuple{String,Any,Symbol}[]
if isfile(ALPHA_MATLAB_COEFFS)
    S = Matrix{Float64}(CSV.read(ALPHA_MATLAB_COEFFS, DataFrame; header=false))
    push!(cands, ("MATLAB coefficients, MATLAB filtfilt", S, :matlab))
    push!(cands, ("MATLAB coefficients, DSP.jl filtfilt",
                  SecondOrderSections([Biquad(r[1] / r[4], r[2] / r[4], r[3] / r[4], r[5] / r[4], r[6] / r[4])
                                       for r in eachrow(S)], 1.0), :zerophase))
end
push!(cands, ("elliptic 6, MATLAB filtfilt", sos_matrix(iir(Elliptic(6, 0.1, 60))), :matlab))
for n in (1, 2, 3, 4, 6, 8)
    push!(cands, ("butterworth $n, zero-phase", iir(Butterworth(n)), :zerophase))
end
push!(cands, ("butterworth 15 (buttord), zero-phase", iir(Butterworth(15), 7.7891, 12.3221), :zerophase))
for n in (4, 6, 8)
    push!(cands, ("chebyshev1 $n, zero-phase", iir(Chebyshev1(n, 0.1)), :zerophase))
end
push!(cands, ("chebyshev2 8, zero-phase", iir(Chebyshev2(8, 60), 6.8497, 13.9905), :zerophase))
for n in (4, 5, 6)
    push!(cands, ("elliptic $n, zero-phase", iir(Elliptic(n, 0.1, 60)), :zerophase))
end
for (name, f) in (("butterworth 2", iir(Butterworth(2))), ("chebyshev1 8", iir(Chebyshev1(8, 0.1))),
                  ("elliptic 6", iir(Elliptic(6, 0.1, 60))))
    push!(cands, ("$name, one pass", f, :causal))
end
for tw in (1.2, 2.0, 4.0)        # FIR, Kaiser window, transition width (Hz)
    push!(cands, ("FIR Kaiser, $(tw) Hz transitions",
                  kaiser_bandpass(fs, 8 - tw / 2, 12 + tw / 2; transition=tw), :fir))
end
push!(cands, ("FIR Kaiser, MATLAB stopbands (6.8, 29.4 Hz)",
              kaiser_bandpass(fs, 7.4, 20.7; transition=1.2), :fir))

# --- Compare ----------------------------------------------------------------------
rows = DataFrame()
for (name, f, mode) in cands
    win_m, win_j, ratio = Float64[], Float64[], Float64[]
    gm = Dict(g => zeros(length(mtimes)) for g in mgroups); gj = deepcopy(gm)
    nsub = Dict(g => 0 for g in mgroups)
    for ((s, g), mp) in mref
        haskey(R, s) || continue
        r = R[s]
        sel = coalesce.(r.condition .== g, false) .& vec(all(!isnan, r.roi; dims=2))
        any(sel) || continue
        P = band_power(r.roi[sel, :], fs, f, mode; pad_s=ALPHA_PAD_S)
        p = vec(mean(P; dims=1))[idx]
        m, j = pct(mp), pct(p)
        push!(win_m, mean(m[win])); push!(win_j, mean(j[win]))
        push!(ratio, median(p) / median(mp))
        gm[g] .+= m; gj[g] .+= j; nsub[g] += 1
    end
    cm = reduce(vcat, [gm[g] ./ nsub[g] for g in mgroups])
    cj = reduce(vcat, [gj[g] ./ nsub[g] for g in mgroups])
    push!(rows, (filter=name, n=length(win_m), r_timecourse=cor(cm, cj),
                 rmse_timecourse=sqrt(mean((cm .- cj) .^ 2)), r_participants=cor(win_m, win_j),
                 window_matlab=mean(win_m), window_julia=mean(win_j), power_ratio=median(ratio)))
end
sort!(rows, :rmse_timecourse)
mkpath(joinpath(RESULTS, "exploratory"))
CSV.write(joinpath(RESULTS, "exploratory", "alpha_filter_identification.csv"), rows)

nsubj = length(unique(first.(collect(keys(mref)))))
println("\nAlpha band-pass filter: Julia candidates vs. MATLAB bandpass() output")
println("(% change from −200–0 ms; window 300–996.1 ms; $(nsubj) participants × 4 conditions)\n")
@printf("%-46s %8s %8s %8s %8s %8s %7s\n", "filter", "r(time)", "RMSE", "r(subj)",
        "MATLAB", "Julia", "power")
for r in eachrow(rows)
    @printf("%-46s %8.3f %8.2f %8.3f %8.1f %8.1f %7.2f\n", r.filter, r.r_timecourse,
            r.rmse_timecourse, r.r_participants, r.window_matlab, r.window_julia, r.power_ratio)
end
println("\nr(time): correlation of the grand-average % change curves; RMSE in % points;")
println("r(subj): correlation of participants' window means; MATLAB/Julia: mean % change")
println("in 300–996.1 ms; power: Julia/MATLAB raw power (1 = same scale).")
println("Best match: $(rows.filter[1]).")
