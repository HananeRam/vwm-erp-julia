# Exploratory analyses: alpha-power suppression and the negative slow wave
# (NSW), following the MATLAB Alpha (Steps 1–3) and NSW (Steps 1–4) scripts.
# Run after 01:  julia --project=. scripts/04_exploratory_alpha_nsw.jl
#
# Outputs in results/exploratory/:
#   alpha_pct_change_by_participant.csv, alpha_tests.csv (t-tests and BF01),
#   alpha_paper_comparison.txt, alpha_figure.png
#   bayes_alpha_<A>_vs_<B>.csv        paired data for Bayesian t-tests in R
#   nsw_cluster_tests.txt, nsw_clusters.csv, nsw_window_tests.csv (t-tests and BF),
#   nsw_paper_comparison.txt, nsw_figure.png
#   bayes_nsw_<window>_<A>_vs_<B>.csv

include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "..", "src", "VWMPipeline.jl"))
include(joinpath(@__DIR__, "..", "src", "ExploratoryStats.jl"))
using .VWMPipeline, .ExploratoryStats
using CSV, DataFrames, Statistics, Plots, DSP

out = joinpath(RESULTS, "exploratory")
mkpath(out)
inwin(t, w) = (t .>= w[1]) .& (t .<= w[2])
safe(s) = replace(s, "-" => "")
"Column means ignoring NaN (participants without trials in a condition)."
colmean(M) = [mean(filter(!isnan, M[:, j])) for j in axes(M, 2)]

"Long table (subject, condition, time, value) → Dict condition => participants × time."
function to_matrices(df, value)
    subs = sort(unique(df.subject))
    times = sort(unique(df.time))
    n_per = combine(groupby(df, :subject), :time => (t -> length(unique(t))) => :n)
    if length(unique(n_per.n)) > 1     # participants on different time grids
        error("Participants have different time points (different sampling rates?): " *
              join(["$(r.subject): $(r.n)" for r in eachrow(n_per)], ", ") *
              ". Rerun scripts/01_process_subjects.jl with the current version.")
    end
    Dict(cond => begin
        M = fill(NaN, length(subs), length(times))
        for r in eachrow(df[df.condition .== cond, :])
            M[findfirst(==(r.subject), subs), findfirst(==(r.time), times)] = r[value]
        end
        M
    end for cond in CONDITIONS), subs, times
end

function wide_table(subs, values::Dict)
    df = DataFrame(subject=subs)
    for cond in CONDITIONS
        df[!, cond] = values[cond]
    end
    return df
end

function export_pairs(wide, pairs, prefix)
    for (a, b) in pairs
        ok = .!isnan.(wide[!, a]) .& .!isnan.(wide[!, b])
        CSV.write(joinpath(out, "$(prefix)_$(safe(a))_vs_$(safe(b)).csv"),
                  DataFrame(:Subject => wide.subject[ok],
                            Symbol(replace(a, "-" => ".")) => wide[ok, a],
                            Symbol(replace(b, "-" => ".")) => wide[ok, b]))
    end
end

# ===========================================================================
# 1. Alpha power: % change from baseline, mean over 300–996.1 ms
# ===========================================================================
alpha_long = CSV.read(joinpath(RESULTS, "alpha.csv"), DataFrame)
alpha_long.condition = String.(alpha_long.condition)
alpha_long.subject = String.(alpha_long.subject)
A, subs, times = to_matrices(alpha_long, :power)
if length(subs) < 3
    @warn "Only $(length(subs)) participant(s) processed: the paired tests and permutation tests need at least 3. Download more data, rerun script 01, then this script."
    exit(0)
end

pct = Dict(c => mapslices(p -> percent_change(p, times, BASELINE), A[c]; dims=2)
           for c in CONDITIONS)
alpha_mean = Dict(c => vec(mean(pct[c][:, inwin(times, ALPHA_WINDOW)]; dims=2))
                  for c in CONDITIONS)
alpha_wide = wide_table(subs, alpha_mean)
CSV.write(joinpath(out, "alpha_pct_change_by_participant.csv"), alpha_wide)

alpha_tests = paired_tests(alpha_wide, PAIRWISE)
CSV.write(joinpath(out, "alpha_tests.csv"), alpha_tests)
export_pairs(alpha_wide, PAIRWISE, "bayes_alpha")
println("\nAlpha power (% change, ", 1000 .* ALPHA_WINDOW, " ms): paired t-tests, Bonferroni")
show(stdout, MIME("text/plain"), alpha_tests); println()

# Comparison with the MATLAB/R Bayes factors (correct trials).
open(joinpath(out, "alpha_paper_comparison.txt"), "w") do io
    println(io, "Alpha power (% change, 300–996.1 ms): Bayesian paired t-tests, BF01")
    println(io, "(Cauchy prior r = √2/2, as BayesFactor::ttestBF in R)\n")
    println(io, rpad("Comparison", 18), rpad("MATLAB", 9), rpad("Julia", 9), "N")
    for (r, (a, b)) in zip(eachrow(alpha_tests), PAIRWISE)
        paper = get(PAPER_ALPHA_BF01, (a, b), NaN)
        println(io, rpad("$a vs $b", 18), rpad(string(round(paper; digits=2)), 9),
                rpad(string(round(r.bf01; digits=2)), 9), r.n)
    end
    println(io, "\nBF01 > 3: moderate evidence for H0; 1–3: anecdotal.")
end
print(read(joinpath(out, "alpha_paper_comparison.txt"), String))

# Figure: time course (grand average, 30 Hz low-pass for display) and bars.
fs_plot = 1 / (times[2] - times[1])
lp = convert(SecondOrderSections,
             digitalfilter(Lowpass(NSW_FIGURE_LOWPASS / (fs_plot / 2)), Butterworth(4)))
# Same style as the paper's Figure 8 (MATLAB).
pa = plot(xlabel="Time (ms)", ylabel="Alpha power (% change from baseline)",
          title="A: Time-resolved alpha power", titlefontsize=11, legend=:topright,
          xlims=(-200, 1000), ylims=(-80, 80), grid=false, framestyle=:box)
vspan!(pa, 1000 .* collect(ALPHA_WINDOW); color=:gray80, alpha=0.4, label="")
vline!(pa, 1000 .* collect(ALPHA_WINDOW); color=:black, ls=:dot, lw=0.8, label="")
for (i, c) in enumerate(CONDITIONS)
    ga = colmean(pct[c])
    plot!(pa, 1000 .* times, filtfilt(lp, ga); lw=2, color=COLORS[i], label=c)
end
vline!(pa, [0]; color=:black, ls=:dash, label="")
ok_vals(c) = filter(!isnan, alpha_wide[!, c])
m = [mean(ok_vals(c)) for c in CONDITIONS]
se = [std(ok_vals(c)) / sqrt(length(ok_vals(c))) for c in CONDITIONS]
ylo = 5 * floor((minimum(m .- se) - 1) / 5)                 # as the MATLAB figure
pb = bar(CONDITIONS, m; yerror=se, color=COLORS, legend=false, grid=false, framestyle=:box,
         ylabel="Mean alpha power suppression (%)", titlefontsize=11,
         title="B: Mean alpha power\n(300–996.1 ms)", ylims=(ylo, 0))
savefig(plot(pa, pb; layout=grid(1, 2; widths=[0.65, 0.35]), size=(1200, 480),
             left_margin=8Plots.mm, bottom_margin=8Plots.mm, top_margin=8Plots.mm),
        joinpath(out, "alpha_figure.png"))

# ===========================================================================
# 2. Negative slow wave: cluster permutation tests and window extraction
# ===========================================================================
nsw_long = CSV.read(joinpath(RESULTS, "nsw.csv"), DataFrame)
nsw_long.condition = String.(nsw_long.condition)
nsw_long.subject = String.(nsw_long.subject)
N, nsubs, ntimes = to_matrices(nsw_long, :amplitude)

cluster_rows = DataFrame()
open(joinpath(out, "nsw_cluster_tests.txt"), "w") do io
    println(io, "Cluster-based permutation tests, NSW (", NSW_PERMUTATIONS,
            " permutations, alpha = ", NSW_ALPHA, ")\n")
    for (a, b) in PAIRWISE
        res = cluster_permutation(N[a], N[b], ntimes; nperm=NSW_PERMUTATIONS, alpha=NSW_ALPHA)
        println(io, "** $a vs $b ** (N = $(res.n), t threshold = $(round(res.t_threshold; digits=3)), ",
                "cluster threshold = $(round(res.cluster_threshold; digits=2)))")
        sig = res.clusters[res.clusters.significant, :]
        if nrow(sig) == 0
            println(io, "  No significant clusters.")
        else
            for r in eachrow(sig)
                println(io, "  $(round(r.start_ms; digits=1)) to $(round(r.end_ms; digits=1)) ms ",
                        "(mass = $(round(r.mass; digits=2)), p = $(round(r.p; digits=4)))")
            end
        end
        println(io)
        res.clusters[!, :comparison] = fill("$a vs $b", nrow(res.clusters))
        append!(cluster_rows, res.clusters; cols=:union)
    end
end
CSV.write(joinpath(out, "nsw_clusters.csv"), cluster_rows)
print(read(joinpath(out, "nsw_cluster_tests.txt"), String))

# Window extraction for the targeted comparisons (Step 4).
window_tests = DataFrame()
for (win, a, b) in NSW_TARGETS
    w = NSW_WINDOWS[win]
    vals = Dict(c => vec(mean(N[c][:, inwin(ntimes, w)]; dims=2)) for c in (a, b))
    wide = DataFrame(:subject => nsubs, Symbol(a) => vals[a], Symbol(b) => vals[b])
    t = paired_tests(wide, [(a, b)]; bonferroni=false)
    t[!, :window] = [win]
    t[!, :window_ms] = ["$(round(1000 * w[1]; digits=1))–$(round(1000 * w[2]; digits=1))"]
    append!(window_tests, t; cols=:union)
    export_pairs(wide, [(a, b)], "bayes_nsw_$(win)")
end
CSV.write(joinpath(out, "nsw_window_tests.csv"), window_tests)
println("\nNSW window tests (uncorrected):")
show(stdout, MIME("text/plain"), window_tests); println()

# Comparison with the MATLAB/R NSW results (correct trials).
fmt_ms(c) = isempty(c) ? "none" : join(["$(round(x; digits=1))–$(round(y; digits=1))" for (x, y) in c], ", ")
open(joinpath(out, "nsw_paper_comparison.txt"), "w") do io
    println(io, "Negative slow wave: Julia vs. MATLAB/R results (correct trials)\n")
    println(io, "Significant clusters (ms)")
    for (a, b) in PAIRWISE
        julia = nrow(cluster_rows) == 0 ? Tuple{Float64,Float64}[] :
                [(r.start_ms, r.end_ms) for r in eachrow(cluster_rows)
                 if r.comparison == "$a vs $b" && r.significant]
        println(io, "  ", rpad("$a vs $b", 16), "MATLAB: ", rpad(fmt_ms(get(PAPER_NSW_CLUSTERS, (a, b), [])), 30),
                "Julia: ", fmt_ms(julia))
    end
    println(io, "\nWindow tests: Bayesian paired t-tests (r = √2/2, as ttestBF in R)")
    println(io, "  ", rpad("Window", 22), rpad("Comparison", 16), rpad("MATLAB BF10", 12), rpad("Julia BF10", 12), "Julia BF01")
    for (r, key) in zip(eachrow(window_tests), NSW_TARGETS)
        fmt(x) = x >= 100 ? string(round(Int, x)) : string(round(x; digits=3))
        println(io, "  ", rpad(key[1], 22), rpad("$(key[2]) vs $(key[3])", 16),
                rpad(fmt(get(PAPER_NSW_BF10, key, NaN)), 12), rpad(fmt(r.bf10), 12), round(r.bf01; digits=3))
    end
    println(io, "\nClusters before 300 ms reflect sensory differences between two- and three-item")
    println(io, "displays; the MATLAB analysis found them too but did not report them.")
    println(io, "\nWindows are the MATLAB cluster windows, so the Julia window tests use the")
    println(io, "same intervals even where the Julia clusters differ slightly.")
end
print(read(joinpath(out, "nsw_paper_comparison.txt"), String))

# Figure: grand average (30 Hz low-pass, display only) with cluster windows.
# Same style as the paper's Figure 7 (MATLAB).
pn = plot(xlabel="Time (ms)", ylabel="ERP amplitude (µV)", titlefontsize=11,
          title="Grand-average negative slow wave (posterior ROI)", legend=:bottomleft,
          grid=false, framestyle=:box)
vspan!(pn, 1000 .* collect(NSW_WINDOWS["Grouping_Cluster"]); color=:gray80, alpha=0.5, label="")
vspan!(pn, 1000 .* collect(NSW_WINDOWS["Late_Cluster"]); color=:gray70, alpha=0.5, label="")
for (i, c) in enumerate(CONDITIONS)
    ga = colmean(N[c])
    # As Step1_Compute_Negative_Slow_Wave.m (figure only): 30 Hz low-pass with
    # pop_eegfiltnew on each epoch, then baseline correction, then averaging.
    # Both steps are linear, so they are applied to each participant's average.
    bl = (ntimes .>= BASELINE[1]) .& (ntimes .<= BASELINE[2])
    rowsf = [begin
                 y = eegfiltnew_lowpass(N[c][k, :], fs_plot, NSW_FIGURE_LOWPASS)
                 y .- mean(y[bl])
             end for k in axes(N[c], 1) if all(!isnan, N[c][k, :])]
    ga = vec(mean(reduce(hcat, rowsf); dims=2))
    plot!(pn, 1000 .* ntimes, ga; lw=2, color=COLORS[i], label=c)
end
vline!(pn, [0]; color=:black, ls=:dash, label="")
plot!(pn; title="A: Grand average negative slow wave", xlims=(-200, 1000), ylims=(-3, 2))

# Panels B and C (as Step3_Plots_Negative_Slow_Wave.m): mean amplitude of the
# unfiltered waves in the two cluster windows, mean ± SE across participants.
function window_bars(win, conds, title)
    vals = [filter(!isnan, vec(mean(N[c][:, inwin(ntimes, NSW_WINDOWS[win])]; dims=2))) for c in conds]
    m = mean.(vals)
    se = std.(vals) ./ sqrt.(length.(vals))
    w = NSW_WINDOWS[win]
    p = bar(conds, m; yerror=se, color=COLORS[indexin(conds, CONDITIONS)], legend=false,
            ylabel="Mean NSW amplitude (µV)", ylims=(-2.5, 1.5), grid=false, framestyle=:box,
            titlefontsize=11,
            title="$title\n($(round(1000w[1]; digits=1))–$(round(1000w[2]; digits=1)) ms)")
    hline!(p, [0]; color=:black, ls=:dash, label="")
    p
end
pB = window_bars("Grouping_Cluster", CONDITIONS, "B: Maintenance cluster")
pC = window_bars("Late_Cluster", ["2-UR", "3-UR"], "C: Late cluster")
savefig(plot(pn, pB, pC; layout=grid(1, 3; widths=[46, 19, 13] ./ 78), size=(1500, 530),
             left_margin=8Plots.mm, bottom_margin=8Plots.mm, top_margin=8Plots.mm),
        joinpath(out, "nsw_figure.png"))

@info "Exploratory results saved to $out"
