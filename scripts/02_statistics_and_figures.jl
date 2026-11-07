# Statistics and figures. Run after 01:
#   julia --project=. scripts/02_statistics_and_figures.jl
#
# 1. Participant level, as in the paper: mean amplitude of each participant's
#    averaged contra–ipsi wave per condition (N2pc, early CDA, sustained CDA);
#    one-way repeated-measures ANOVAs and pairwise comparisons; a table
#    comparing the results with the published values.
# 2. Single-trial level (new): linear mixed models with participant random effects.

include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "..", "src", "ExploratoryStats.jl"))
using .ExploratoryStats
include(joinpath(@__DIR__, "..", "src", "VWMPipeline.jl"))
using .VWMPipeline: erplab_window
using CSV, DataFrames, Statistics, MixedModels, Effects, Plots, Printf

trials = CSV.read(joinpath(RESULTS, "trials.csv"), DataFrame)
waves  = CSV.read(joinpath(RESULTS, "waveforms.csv"), DataFrame)
heog   = CSV.read(joinpath(RESULTS, "heog.csv"), DataFrame)
qc     = CSV.read(joinpath(RESULTS, "qc.csv"), DataFrame)
for df in (trials, waves, heog)
    df.condition = String.(df.condition)
    df.subject   = String.(df.subject)
end
kept = trials[trials.keep, :]
contrasts = Dict(:condition => DummyCoding(base="2-UR"))

# --- 1. Participant-level measurements (Script #11a, pop_geterpvalues) -------
# Mean amplitude ('meanbl', baseline −200–0 ms) on the averaged waves WITHOUT
# the 20 Hz low-pass; window edges at the closest samples, as in ERPLAB.
sort!(waves, [:subject, :condition, :time])
rebase(t, a) = a .- mean(a[(t .>= BASELINE[1]) .& (t .<= BASELINE[2])])
meanbl(t, a, w) = mean(rebase(t, a)[erplab_window(t, w)])
subj = combine(groupby(waves, [:subject, :condition]),
    [:time, :amplitude_raw] => ((t, a) -> meanbl(t, a, N2PC_WINDOW)) => :n2pc,
    [:time, :amplitude_raw] => ((t, a) -> meanbl(t, a, CDA_EARLY_WINDOW)) => :cda_early,
    [:time, :amplitude_raw] => ((t, a) -> meanbl(t, a, CDA_WINDOW))  => :cda,
    :n_trials => first => :n_trials)
CSV.write(joinpath(RESULTS, "participant_means.csv"), subj)

measures = (:n2pc, :cda_early, :cda)
cond_summary = combine(groupby(subj, :condition), nrow => :n,
    [m => mean => Symbol(m, "_mean") for m in measures]...,
    [m => std => Symbol(m, "_sd") for m in measures]...,
    [m => (x -> std(x) / sqrt(length(x))) => Symbol(m, "_se") for m in measures]...)
cond_summary = cond_summary[indexin(CONDITIONS, cond_summary.condition), :]
CSV.write(joinpath(RESULTS, "condition_summary.csv"), cond_summary)

# --- Repeated-measures ANOVAs and pairwise comparisons (as in the paper) ------
wide_of(m) = unstack(subj, :subject, :condition, m)
anova_rows = DataFrame()
pairwise_all = DataFrame()      # pooled SE, as in the R analysis (lme + glht)
paired_all = DataFrame()        # paired t-tests, for reference
enough = length(unique(subj.subject)) >= 3      # ANOVA needs several participants
enough || @warn "ANOVAs skipped: fewer than 3 participants."
for m in (enough ? measures : ())
    w = wide_of(m)
    a = rm_anova(w, CONDITIONS)
    push!(anova_rows, (measure=String(m), n=a.n, F=a.F, df1=a.df1, df2=a.df2,
                       p=a.p, eta_p2=a.eta_p2); cols=:union)
    pw = pairwise_pooled(w, ERP_PAIRS, a.ms_error, a.df2)
    pw[!, :measure] = fill(String(m), nrow(pw))
    append!(pairwise_all, pw; cols=:union)
    pt = pairwise_comparisons(w, ERP_PAIRS)
    pt[!, :measure] = fill(String(m), nrow(pt))
    append!(paired_all, pt; cols=:union)
end
CSV.write(joinpath(RESULTS, "erp_anova.csv"), anova_rows)
CSV.write(joinpath(RESULTS, "erp_pairwise.csv"), pairwise_all)
CSV.write(joinpath(RESULTS, "erp_pairwise_paired_t.csv"), paired_all)

# --- Comparison with the published results ----------------------------------
if enough
open(joinpath(RESULTS, "paper_comparison.txt"), "w") do io
    labels = Dict("n2pc" => "N2pc (238.3–300 ms)", "cda_early" => "Early CDA (300–450 ms)",
                  "cda" => "Sustained CDA (300–996.1 ms)")
    println(io, "Julia reproduction vs. MATLAB/R results (correct trials, participant-level means)\n")
    for m in ("n2pc", "cda_early", "cda")
        println(io, labels[m])
        a = anova_rows[anova_rows.measure .== m, :]
        paper = get(PAPER_ANOVA, m, nothing)
        paper === nothing ||
            @printf(io, "  ANOVA  MATLAB: F(%d, %d) = %.2f, p %s, ηp² = %.3f\n",
                    paper.df1, paper.df2, paper.F, paper.p, paper.eta_p2)
        @printf(io, "         Julia:  F(%d, %d) = %.2f, p = %.4f, ηp² = %.3f (N = %d)\n",
                a.df1[1], a.df2[1], a.F[1], a.p[1], a.eta_p2[1], a.n[1])
        println(io, "  Mean (SD), µV      MATLAB         Julia")
        for c in CONDITIONS
            r = cond_summary[cond_summary.condition .== c, :]
            jm, jsd = r[1, Symbol(m, "_mean")], r[1, Symbol(m, "_sd")]
            pv = haskey(PAPER_MEANS, m) ? PAPER_MEANS[m][c] : nothing
            ptxt = pv === nothing ? "      —      " : @sprintf("%6.2f (%4.2f)", pv[1], pv[2])
            @printf(io, "    %-7s      %s   %6.2f (%4.2f)\n", c, ptxt, jm, jsd)
        end
        println(io, "  Pairwise (A − B), Bonferroni   MATLAB: diff [95% CI], p  Julia: diff [95% CI], p, dz")
        for r in eachrow(pairwise_all[pairwise_all.measure .== m, :])
            a_, b_ = strip.(split(r.comparison, "−"))
            pp = get(get(PAPER_PAIRWISE, m, Dict()), (a_, b_), nothing)
            ptxt = pp === nothing ? rpad("not reported", 28) :
                   isnan(pp[1]) ? rpad(@sprintf("  —   p = %.3f", pp[4]), 28) :
                   @sprintf("%5.2f [%5.2f, %5.2f], p = %.3f", pp[1], pp[2], pp[3], pp[4])
            @printf(io, "    %-15s %s   %5.2f [%5.2f, %5.2f], p = %.3f, dz = %.2f\n",
                    r.comparison, ptxt, r.mean_diff, r.ci_low, r.ci_high,
                    r.p_bonferroni, r.cohens_dz)
        end
        println(io)
    end
    println(io, "Pairwise method as in the R analysis (lme + glht, Bonferroni over six comparisons):")
    println(io, "pooled SE = √(2·MS_error/n), z-tests, 95% CI = diff ± t(.975, df_error)·SE.")
    println(io, "Uncorrected p in erp_pairwise.csv; paired t-tests in erp_pairwise_paired_t.csv.")
    println(io, "MATLAB values: ERPLAB + R analysis (correct trials only), as reported")
    println(io, "in the Stage 2 manuscript (2 decimals); \"—\" / not reported = not in the text.")
end
print(read(joinpath(RESULTS, "paper_comparison.txt"), String))
end

# --- 2. Mixed models -----------------------------------------------------------
# Mixed models need several participants: with fewer than 3, a participant
# random effect cannot be estimated, so the models are skipped.
results = Dict{String,Any}()
n_subjects = length(unique(subj.subject))
open(joinpath(RESULTS, "models.txt"), "w") do io
    if n_subjects < 3
        msg = "Mixed models skipped: only $n_subjects participant(s). Download more data."
        println(io, msg)
        @warn msg
        return
    end
    for (label, df) in (("participant means", subj), ("single trials", kept)),
        dv in (:cda, :n2pc)
        f = dv === :cda ? @formula(cda ~ 1 + condition + (1 | subject)) :
                          @formula(n2pc ~ 1 + condition + (1 | subject))
        println(io, "=== ", uppercase(String(dv)), " (", label, ") ===")
        try
            m = fit(MixedModel, f, df; contrasts=contrasts, progress=false)
            results["$(dv) $(label)"] = m
            show(io, MIME("text/plain"), m)
        catch err
            println(io, "Model failed: ", sprint(showerror, err))
            @warn "Model for $(dv) ($(label)) failed" exception=err
        end
        println(io, "\n")
    end
end

if haskey(results, "cda single trials")
    eff = effects(Dict(:condition => CONDITIONS), results["cda single trials"])
    CSV.write(joinpath(RESULTS, "cda_effects.csv"), eff)
    p1 = scatter(eff.condition, eff.cda; yerror=eff.err, legend=false,
                 xlabel="Condition", ylabel="CDA amplitude (µV)",
                 title="CDA by condition (single-trial mixed model ± SE)")
    savefig(p1, joinpath(RESULTS, "cda_effects.png"))
end

# --- Figures -------------------------------------------------------------------
# Same style as the paper's Figures 3 and 4 (MATLAB): ms axis, the paper's colors,
# window arrows below the waves, electrode names top left.
function grand_average_plot(df, title, ylabel; ylims=(-1.2, 1), zero_line_top=0.85)
    ga = combine(groupby(df, [:condition, :time]), :amplitude => mean => :amplitude)
    p = plot(xlabel="Time from array onset (ms)", ylabel=ylabel, title=title,
             xlims=(-200, 1000), ylims=ylims, legend=:topright, foreground_color_legend=nothing,
             grid=false, size=(1000, 420), left_margin=5Plots.mm, bottom_margin=5Plots.mm)
    hline!(p, [0]; color=:black, lw=0.8, label="")
    # Dashed 0-ms line stops below the electrode names so it does not cross them.
    plot!(p, [0, 0], [ylims[1], zero_line_top]; color=:black, ls=:dash, lw=0.8, label="")
    for (i, cond) in enumerate(CONDITIONS)
        d = ga[ga.condition .== cond, :]
        plot!(p, 1000 .* d.time, d.amplitude; lw=2, color=COLORS[i], label=cond)
    end
    return p
end

# Window arrow with its label underneath (as in the MATLAB figures).
function window_arrow!(p, win_s, label; y=-0.93, gap_ms=4)
    a, b = 1000 .* collect(win_s) .+ [gap_ms, -gap_ms] ./ 2
    plot!(p, [a, b], [y, y]; color=:black, lw=1.5, label="")
    scatter!(p, [a], [y]; marker=:ltriangle, ms=4, color=:black, markerstrokewidth=0, label="")
    scatter!(p, [b], [y]; marker=:rtriangle, ms=4, color=:black, markerstrokewidth=0, label="")
    annotate!(p, (a + b) / 2, y - 0.04, text(label, 10, :center, :top))
    return p
end

const ELECTRODE_LABEL = "P3/4, P5/6, P7/8, P9/10, PO3/4, PO7/8, O1/2"

# Figure 4 style: N2pc and sustained CDA windows.
p2 = grand_average_plot(waves, "", "Contra − ipsi amplitude (µV)"; ylims=(-1.2, 1))
window_arrow!(p2, N2PC_WINDOW, "N2pc")
window_arrow!(p2, CDA_WINDOW, "Sustained CDA")
annotate!(p2, -185, 0.97, text(ELECTRODE_LABEL, 10, :left, :top))
savefig(p2, joinpath(RESULTS, "grand_average.png"))

# Figure 3 style: N2pc and early CDA windows, −200 to 500 ms.
p2b = grand_average_plot(waves, "", "Contra − ipsi amplitude (µV)"; ylims=(-1.2, 1))
plot!(p2b; xlims=(-200, 500), size=(800, 420))
window_arrow!(p2b, N2PC_WINDOW, "N2pc")
window_arrow!(p2b, CDA_EARLY_WINDOW, "Early CDA")
annotate!(p2b, -190, 0.97, text(ELECTRODE_LABEL, 10, :left, :top))
savefig(p2b, joinpath(RESULTS, "grand_average_early.png"))

# Eye-movement check (Script #8a): residual HEOG should be small (Luck: < 3.2 µV).
p3 = grand_average_plot(heog, "HEOG contra − ipsi (eye-movement check)", "HEOG contra − ipsi (µV)";
                        ylims=(-1, 1), zero_line_top=1)
savefig(p3, joinpath(RESULTS, "heog_check.png"))

p4 = bar(1:nrow(qc), qc.pct_rejected_total; legend=false, color=COLORS[1], grid=false,
         xticks=(1:nrow(qc), replace.(qc.subject, "sub-" => "")), xrotation=90, xtickfontsize=7,
         xlabel="Participant", ylabel="% trials rejected", title="Rejected trials per participant",
         size=(900, 400), left_margin=5Plots.mm, bottom_margin=5Plots.mm)
savefig(p4, joinpath(RESULTS, "rejection_rates.png"))

@info "Saved tables, models.txt and figures to $RESULTS"
