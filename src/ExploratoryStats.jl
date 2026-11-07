"""
    ExploratoryStats

Paired tests used in the exploratory alpha-power and negative-slow-wave
analyses, written to match the MATLAB scripts:

- `paired_tests`: paired t-tests with Bonferroni correction and Cohen's dz
  (Step3_Bays_Factor_Alpha_Power, Step4_Bayesian_Window_Extraction_NSW)
- `jzs_bf10`: Bayes factor for a paired t-test with a Cauchy prior on effect
  size (r = √2/2), as `BayesFactor::ttestBF(x, y, paired = TRUE)` in R
- `cluster_permutation`: cluster-based permutation test on paired difference
  waves with sign flipping, cluster mass = sum of |t|, max-cluster null
  distribution (Step2_Permutation_Test_Negative_Slow_Wave)
"""
module ExploratoryStats

using Statistics, Random, Distributions, HypothesisTests, DataFrames

export jzs_bf10, paired_tests, cluster_permutation, contiguous_runs, rm_anova, pairwise_comparisons,
       pairwise_pooled

"""
    rm_anova(wide, conditions) -> NamedTuple

One-way repeated-measures ANOVA on a participants × conditions table
(complete cases only): F, df, p, and partial eta squared
(SS_condition / (SS_condition + SS_error)). No sphericity correction.
"""
function rm_anova(wide::DataFrame, conditions)
    X = Matrix{Float64}(coalesce.(wide[:, conditions], NaN))
    X = X[vec(all(!isnan, X; dims=2)), :]
    n, k = size(X)
    grand = mean(X)
    ss_cond = n * sum((vec(mean(X; dims=1)) .- grand) .^ 2)
    ss_subj = k * sum((vec(mean(X; dims=2)) .- grand) .^ 2)
    ss_total = sum((X .- grand) .^ 2)
    ss_error = ss_total - ss_cond - ss_subj
    df1, df2 = k - 1, (k - 1) * (n - 1)
    F = (ss_cond / df1) / (ss_error / df2)
    return (n=n, F=F, df1=df1, df2=df2, p=ccdf(FDist(df1, df2), F),
            eta_p2=ss_cond / (ss_cond + ss_error), ms_error=ss_error / df2)
end

"""
    pairwise_pooled(wide, pairs, ms_error, df_error) -> DataFrame

Pairwise comparisons as in the original R analysis: `lme(y ~ condition,
random = ~1 | subject/condition)` followed by `glht(..., mcp(condition =
"Tukey"))` with Bonferroni adjustment. With one value per participant and
condition, the standard error of every difference is the same,
√(2 · MS_error / n), taken from the repeated-measures ANOVA. Tests are z-tests
(as in `glht`); 95% CIs use t(0.975, df_error), as in the R script.
"""
function pairwise_pooled(wide::DataFrame, pairs, ms_error::Real, df_error::Integer)
    conds = unique(vcat(collect.(pairs)...))
    X = coalesce.(wide[:, conds], NaN)
    complete = [all(c -> !isnan(X[i, c]), conds) for i in 1:nrow(X)]   # same rows as rm_anova
    n = count(complete)
    se = sqrt(2 * ms_error / n)
    tc = quantile(TDist(df_error), 0.975)
    rows = map(pairs) do (a, b)
        d = Float64.(X[complete, a]) .- Float64.(X[complete, b])
        md = mean(d)
        z = md / se
        p = 2 * ccdf(Normal(), abs(z))
        (comparison="$a − $b", n=n, mean_diff=md, se=se, ci_low=md - tc * se,
         ci_high=md + tc * se, z=z, p=p, p_bonferroni=min(1.0, p * length(pairs)),
         cohens_dz=md / std(d))
    end
    return DataFrame(rows)
end

"""
    jzs_bf10(t, n; r=√2/2) -> Float64

JZS Bayes factor BF10 for a one-sample or paired t-test (Rouder et al., 2009),
with a Cauchy prior of scale `r` on the effect size; the default matches the
"medium" prior of `BayesFactor::ttestBF` in R. BF01 = 1 / BF10.
The integral over g is computed on a log scale (trapezoid rule).
Check: R's sleep data (paired, n = 10, t = −4.06) gives BF10 = 17.2589, the
value in the BayesFactor documentation.
"""
function jzs_bf10(t::Real, n::Integer; r::Real=sqrt(2) / 2)
    ν = n - 1
    us = range(-30.0, 30.0; length=60_001)          # u = log(g)
    integrand(u) = begin
        g = exp(u)
        a = 1 + n * r^2 * g
        a^-0.5 * (1 + t^2 / (a * ν))^(-(ν + 1) / 2) *
            (2π)^-0.5 * g^-1.5 * exp(-1 / (2g)) * g  # × g: dg = g du
    end
    v = integrand.(us)
    integral = step(us) * (sum(v) - (first(v) + last(v)) / 2)
    return integral / (1 + t^2 / ν)^(-(ν + 1) / 2)
end

"Holm step-down adjustment of a vector of p-values."
function holm(p::AbstractVector)
    m = length(p)
    order = sortperm(p)
    adj = similar(p, Float64)
    running = 0.0
    for (rank, i) in enumerate(order)
        running = max(running, min(1.0, (m - rank + 1) * p[i]))
        adj[i] = running
    end
    return adj
end

"""
    pairwise_comparisons(wide, pairs) -> DataFrame

Paired comparisons A − B (complete pairs): mean difference with 95% CI,
t, df, uncorrected p, Bonferroni and Holm adjusted p (over all `pairs`),
and Cohen's dz (mean difference / SD of differences).
"""
function pairwise_comparisons(wide::DataFrame, pairs)
    rows = map(pairs) do (a, b)
        x, y = coalesce.(wide[!, a], NaN), coalesce.(wide[!, b], NaN)
        ok = .!isnan.(x) .& .!isnan.(y)
        d = Float64.(x[ok] .- y[ok])
        n = length(d)
        se = std(d) / sqrt(n)
        tc = quantile(TDist(n - 1), 0.975)
        t = mean(d) / se
        (comparison="$a − $b", n=n, mean_diff=mean(d), ci_low=mean(d) - tc * se,
         ci_high=mean(d) + tc * se, t=t, df=n - 1, p=2 * ccdf(TDist(n - 1), abs(t)),
         cohens_dz=mean(d) / std(d))
    end
    df = DataFrame(rows)
    df.p_bonferroni = min.(1.0, df.p .* nrow(df))
    df.p_holm = holm(df.p)
    return df
end

"""
    paired_tests(wide, pairs; bonferroni=true) -> DataFrame

`wide` has one row per participant and one column per condition. For each
pair (A, B): paired t-test on participants with both values, t, df, p,
Bonferroni-corrected p (p × number of pairs, capped at 1), Cohen's dz, and
Bayes factors BF10 and BF01 (`jzs_bf10`, as `ttestBF` in R).
"""
function paired_tests(wide::DataFrame, pairs; bonferroni=true)
    m = length(pairs)
    rows = map(pairs) do (a, b)
        ok = .!ismissing.(wide[!, a]) .& .!ismissing.(wide[!, b]) .&
             .!isnan.(coalesce.(wide[!, a], NaN)) .& .!isnan.(coalesce.(wide[!, b], NaN))
        x, y = Float64.(wide[ok, a]), Float64.(wide[ok, b])
        if length(x) < 3                               # too few participants to test
            return (comparison="$a vs $b", n=length(x), mean_diff=NaN, t=NaN, df=NaN,
                    p=NaN, p_bonferroni=NaN, cohens_dz=NaN, bf10=NaN, bf01=NaN)
        end
        t = OneSampleTTest(x, y)                       # paired t-test
        d = x .- y
        p = pvalue(t)
        bf = jzs_bf10(t.t, length(x))
        (comparison="$a vs $b", n=length(x), mean_diff=mean(d), t=t.t, df=t.df,
         p=p, p_bonferroni=bonferroni ? min(1.0, p * m) : p, cohens_dz=mean(d) / std(d),
         bf10=bf, bf01=1 / bf)
    end
    return DataFrame(rows)
end

"Index ranges of consecutive `true` values in a Boolean vector."
function contiguous_runs(mask::AbstractVector{Bool})
    runs = UnitRange{Int}[]
    i, n = 1, length(mask)
    while i <= n
        if mask[i]
            j = i
            while j < n && mask[j+1]
                j += 1
            end
            push!(runs, i:j)
            i = j + 1
        else
            i += 1
        end
    end
    return runs
end

paired_t(D::AbstractMatrix) =
    vec(mean(D; dims=1) ./ (std(D; dims=1) ./ sqrt(size(D, 1))))

"""
    cluster_permutation(A, B, times; nperm, alpha, rng) -> NamedTuple

`A` and `B` are participants × time matrices for two conditions. Points with
|t| above the two-tailed critical t form clusters; cluster mass is the sum of
|t|. The null distribution is the maximum cluster mass after randomly
flipping the sign of each participant's difference wave. A cluster is
significant if its mass exceeds the (1 − alpha) percentile of the null.
"""
function cluster_permutation(A::AbstractMatrix, B::AbstractMatrix, times::AbstractVector;
                             nperm::Integer=10_000, alpha::Real=0.05,
                             rng::AbstractRNG=MersenneTwister(2026))
    valid = [all(!isnan, A[i, :]) && all(!isnan, B[i, :]) for i in axes(A, 1)]
    D = A[valid, :] .- B[valid, :]
    n = size(D, 1)
    tcrit = quantile(TDist(n - 1), 1 - alpha / 2)

    tobs = paired_t(D)
    obs = [(r, sum(abs, tobs[r])) for r in contiguous_runs(abs.(tobs) .> tcrit)]

    null = zeros(nperm)
    for p in 1:nperm
        signs = rand(rng, (-1.0, 1.0), n)
        tp = paired_t(D .* signs)
        rr = contiguous_runs(abs.(tp) .> tcrit)
        null[p] = isempty(rr) ? 0.0 : maximum(sum(abs, tp[r]) for r in rr)
    end
    cthr = quantile(null, 1 - alpha)

    clusters = DataFrame(start_ms=[1000 * times[first(r)] for (r, _) in obs],
                         end_ms=[1000 * times[last(r)] for (r, _) in obs],
                         mass=[m for (_, m) in obs],
                         p=[mean(null .>= m) for (_, m) in obs])
    clusters.significant = clusters.mass .> cthr
    return (n=n, t_threshold=tcrit, cluster_threshold=cthr, clusters=clusters, t=tobs)
end

end # module
