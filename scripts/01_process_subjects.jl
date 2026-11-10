# Process every listed participant whose files are present. Mirrors MATLAB
# Scripts #1, #5, #6a and #7a (without ICA; see README). Saves:
#   results/trials.csv      one row per trial: labels, artifact flags, N2pc, CDA
#   results/waveforms.csv   participant × condition cluster contra–ipsi waves
#   results/heog.csv        participant × condition HEOG contra–ipsi waves
#   results/qc.csv          trials and % rejected per participant and detector
#
# Run:  julia --project=. scripts/01_process_subjects.jl

include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "..", "src", "VWMPipeline.jl"))
using .VWMPipeline
using CSV, DataFrames, Statistics, DSP, Serialization

eeg_path(sub)    = joinpath(BIDS_ROOT, sub, "eeg", "$(sub)_task-$(TASK)_eeg.vhdr")
bin_path(sub)    = joinpath(BIDS_ROOT, sub, "eeg", "$(sub)_task-$(TASK)_eeg.eeg")
events_path(sub) = joinpath(BIDS_ROOT, sub, "eeg", "$(sub)_task-$(TASK)_events.tsv")

# The download loop fetches events.tsv last, so its presence means the
# participant's files are complete.
downloaded(d) = isfile(eeg_path(d)) && isfile(bin_path(d)) &&
                filesize(bin_path(d)) > 1_000_000 && isfile(events_path(d))
subjects = filter(downloaded, SUBJECTS)
isempty(subjects) && error("No EEG files found under $BIDS_ROOT. See README: Get the data.")
missing_subs = setdiff(SUBJECTS, subjects)
isempty(missing_subs) ||
    @warn "$(length(missing_subs)) of $(length(SUBJECTS)) listed participants not downloaded yet: $(join(missing_subs, ", "))"

# Optional per-participant thresholds (id, name, threshold).
overrides = isfile(ARTIFACT_OVERRIDES) ? CSV.read(ARTIFACT_OVERRIDES, DataFrame) : nothing
function threshold_for(id, art)
    overrides === nothing && return art.threshold
    r = overrides[(overrides.id .== id) .& (overrides.name .== art.name), :]
    return nrow(r) == 0 ? art.threshold : Float64(r.threshold[1])
end
const DETECTORS = Dict(:step => art_step, :extreme => art_extreme, :mwpp => art_mwpp)
cluster_channels = vcat([[l, r] for (l, r) in CLUSTER]...)

mkpath(RESULTS)
@info "Blink correction: $(BLINK_CORRECTION); alpha band-pass filter: $(ALPHA_FILTER)"
trials_all, waves_all, heog_all, qc_all, nsw_all, alpha_all, ica_all =
    (DataFrame() for _ in 1:7)
alpha_cand = Dict(d => DataFrame() for d in ALPHA_CANDIDATES)
roi_trials = Dict{String,Any}()          # ROI trials for the alpha filter search (script 05)

for sub in subjects
    id = parse(Int, sub[5:end])
    @info "Processing $sub"

    # Script #1: load, downsample to 256 Hz, re-reference, bipolar EOG, DC
    # removal, high-pass.
    data, fs = load_eeg(eeg_path(sub); channels=USE_ICA ? nothing : LOAD_CHANNELS)
    fs_orig = fs
    if fs != TARGET_FS                 # e.g. a recording saved at 1024 Hz
        r = rationalize(TARGET_FS / fs)
        for ch in collect(keys(data))
            # Script #1, pop_resample: EEGLAB's own anti-aliasing FIR (integer factors);
            # DSP.jl's resample otherwise.
            data[ch] = numerator(r) == 1 ? eeglab_resample(data[ch], denominator(r)) :
                                           resample(Float64.(data[ch]), r)
        end
        fs = TARGET_FS
    end
    if BLINK_CORRECTION === :ica_weights
        # The original ICA was computed on data referenced to the left mastoid
        # (Script #1); the weights are applied in that reference, then the data
        # are re-referenced to the average of the mastoids (Script #4).
        lm = copy(data["Mastoid-left"])
        for (_, x) in data
            x .-= lm
        end
    else
        rereference!(data)
    end
    add_bipolars!(data)
    remove_dc!(data)

    butter_filter!(data, fs; kind=:highpass, cutoff=HIGHPASS.cutoff, order=HIGHPASS.order)

    # Script #5: trial labels and epochs.
    events = read_events(events_path(sub))
    # sub-10 and sub-49: scripts/raw_events/<sub>.csv holds the events of the MATLAB
    # Script #0 dataset (1-based latencies at 512 Hz; 1248 and 1244 memory arrays),
    # the input of the original analysis. Since dataset version 1.1.0 the OpenNeuro
    # files of these two participants are the same recordings and events; in versions
    # 1.0.x the events were incomplete and sub-49 was at 2048 Hz (README).
    # Latencies are rescaled if the EEG file has another rate.
    raw_ev = joinpath(@__DIR__, "raw_events", "$(sub).csv")
    if isfile(raw_ev)
        r = CSV.read(raw_ev, DataFrame; types=Dict(:typ => String))
        val = [something(tryparse(Int, filter(isdigit, t)), -1) for t in r.typ]
        smp = round.(Int, (r.lat .- 1) .* (fs_orig / 512)) .+ 1
        events = DataFrame(sample=smp, value=val)
    end
    if fs_orig != fs                   # event samples (1-based) to the new rate
        events.sample = floor.(Int, (events.sample .- 1) .* (fs / fs_orig)) .+ 1
    end
    # Scripts #2–#4: ICA blink correction of the scalp channels (EOG channels and
    # the bipolar HEOG/VEOG stay uncorrected for artifact detection).
    if BLINK_CORRECTION === :ica_weights
        rm_ics = get(ICS_APPLIED, id, ics_to_remove(ICS_TO_REMOVE, id))
        ica = read_ica_weights(find_ica_set(ICA_WEIGHTS_DIR, id), rm_ics)
        apply_ica_weights!(data, ica)
        rereference!(data)                      # left mastoid → average of mastoids
        push!(ica_all, (subject=sub, n_channels=length(ica.channels),
                        n_components=size(ica.unmixing, 1), n_removed=length(ica.remove),
                        components=join(ica.remove, ";")))
    elseif BLINK_CORRECTION === :fastica
        scalp = sort([ch for ch in keys(data) if !(ch in ("Mastoid-left", "Mastoid-right",
                      "HEOG-left", "HEOG-right", "VEOG-lower", "HEOG", "VEOG"))])
        ica = ica_blink_correct!(data, fs, collect(skipmissing(events.sample));
                                 scalp=scalp, n_remove=ica_n_remove(id),
                                 n_components=ICA_COMPONENTS, seed=id)
        push!(ica_all, (subject=sub, n_channels=length(scalp), n_removed=length(ica.removed),
                        components=join(ica.removed, ";"),
                        r_veog=join(round.(ica.r; digits=2), ";")))
    end
    if USE_ICA
        for ch in setdiff(keys(data), vcat(LOAD_CHANNELS, ["HEOG", "VEOG"]))
            delete!(data, ch)            # free memory: keep only the channels analyzed
        end
        GC.gc()
    end
    codes = event_codes(events)
    id in SWAP_RESPONSE_IDS && swap_response_codes!(codes)
    trials = label_trials(codes, events.sample)
    onsets = trials.stim_sample
    epoch_all(d) = Dict(ch => first(epoch(d[ch], onsets, fs; tmin=EPOCH.tmin,
                                          tmax=EPOCH.tmax, baseline=BASELINE))
                        for ch in vcat(cluster_channels, ["HEOG", "VEOG"]))
    times = last(epoch(data["HEOG"], onsets[1:1], fs; tmin=EPOCH.tmin, tmax=EPOCH.tmax))
    ep_art = epoch_all(data)          # high-pass only: used for artifact detection

    # Script #7a's 20 Hz low-pass on the epochs (= on the averaged ERPs, as in
    # pop_filterp), used for ERP measurement. Artifact detection, the NSW and
    # alpha use the unfiltered epochs (ep_art), as in the MATLAB scripts.
    ep = Dict(ch => lowpass_epochs!(copy(E), fs; cutoff=LOWPASS.cutoff,
                                    order=LOWPASS.order) for (ch, E) in ep_art)

    # Script #6a: artifact detection, one flag column per detector.
    # Rates are relative to binned trials (correct responses), as in ERPLAB.
    binned = CORRECT_ONLY ? coalesce.(trials.correct, false) : trues(nrow(trials))
    any_artifact = falses(nrow(trials))
    qc = Dict{Symbol,Any}(:subject => sub, :n_trials_total => nrow(trials),
                          :n_binned => count(binned),
                          :n_inconsistent_codes => count(coalesce.(trials.consistent, false) .== false))
    for art in ARTIFACTS
        det = DETECTORS[art.detector]
        thr = threshold_for(id, art)
        chans = art.channel === :cluster ? cluster_channels : [art.channel]
        flagged = falses(nrow(trials))
        for ch in chans
            flagged .|= flag_trials(ep_art[ch], times, fs, det; period=art.period,
                                    threshold=thr, window_ms=art.window_ms,
                                    step_ms=art.step_ms)
        end
        trials[!, Symbol("art_", art.name)] = flagged
        qc[Symbol("pct_", art.name)] = round(100 * mean(flagged[binned]); digits=1)
        any_artifact .|= flagged
    end
    trials.artifact = any_artifact
    keep = binned .& .!any_artifact
    trials.keep = keep
    qc[:n_kept] = count(keep)
    qc[:pct_rejected_total] = round(100 * (1 - count(keep) / count(binned)); digits=1)
    push!(qc_all, NamedTuple(qc); cols=:union)

    # Script #7a: cluster contra-minus-ipsi, average waves.
    # Script #11a measures mean amplitudes on the averages WITHOUT the 20 Hz
    # low-pass (D_raw); the low-passed waves (D) are for the figures.
    D = cluster_contra_minus_ipsi(ep, CLUSTER, trials.cue)
    D_raw = cluster_contra_minus_ipsi(ep_art, CLUSTER, trials.cue)
    H = similar(ep["HEOG"])                      # HEOG contra − ipsi (sign by cue)
    for k in axes(H, 1)
        H[k, :] .= trials.cue[k] == "left" ? .-ep["HEOG"][k, :] : ep["HEOG"][k, :]
    end
    trials.n2pc = mean_amplitude(D_raw, times, N2PC_WINDOW)
    trials.cda  = mean_amplitude(D_raw, times, CDA_WINDOW)
    trials.subject = fill(sub, nrow(trials))
    append!(trials_all, select(trials, :subject, :condition, :cue, :change, :probe,
                               :correct, :consistent, r"^art_", :keep, :n2pc, :cda); cols=:union)

    # Exploratory: non-lateralized ROI average (14 channels, no 20 Hz low-pass)
    # for the negative slow wave, and its time-resolved alpha power.
    roi = sum(ep_art[c] for c in cluster_channels) ./ length(cluster_channels)
    roi_kept = fill(NaN, size(roi))
    roi_kept[keep, :] .= roi[keep, :]
    roi_trials[sub] = (times=times, roi=roi_kept, condition=copy(trials.condition), fs=fs)
    alpha_by = Dict(d => alpha_power(roi_kept, fs; band=ALPHA_BAND, pad_s=ALPHA_PAD_S,
                                     hilbert_padded=ALPHA_HILBERT_PADDED, design=d)
                    for d in ALPHA_CANDIDATES)
    alpha = alpha_power(roi_kept, fs; band=ALPHA_BAND, pad_s=ALPHA_PAD_S,
                        hilbert_padded=ALPHA_HILBERT_PADDED, design=ALPHA_FILTER,
                        matlab_coeffs=ALPHA_MATLAB_COEFFS)

    for cond in CONDITIONS
        sel = keep .& (trials.condition .== cond)
        any(sel) || continue
        avg(M) = vec(mean(M[sel, :]; dims=1))
        # Lateralized waves as in the ERPLAB bin operations (Script #7a):
        # left-cue and right-cue trials are averaged separately, then combined
        # 50/50 (nbin = 0.5*bin_left + 0.5*bin_right), so each hemifield has
        # the same weight whatever its number of kept trials.
        sel_l = sel .& (trials.cue .== "left")
        sel_r = sel .& (trials.cue .== "right")
        (any(sel_l) && any(sel_r)) ||
            @warn "Participant $sub, $cond: no kept trials on one side ($(count(sel_l)) left, $(count(sel_r)) right)"
        avg_lat(M) = !any(sel_l) ? vec(mean(M[sel_r, :]; dims=1)) :
                     !any(sel_r) ? vec(mean(M[sel_l, :]; dims=1)) :
                     0.5 .* (vec(mean(M[sel_l, :]; dims=1)) .+ vec(mean(M[sel_r, :]; dims=1)))
        append!(waves_all, DataFrame(subject=sub, condition=cond, n_trials=count(sel),
                                     n_left=count(sel_l), n_right=count(sel_r),
                                     time=times, amplitude=avg_lat(D), amplitude_raw=avg_lat(D_raw)))
        append!(heog_all, DataFrame(subject=sub, condition=cond, time=times, amplitude=avg_lat(H)))
        append!(nsw_all, DataFrame(subject=sub, condition=cond, n_trials=count(sel),
                                   time=times, amplitude=avg(roi)))
        append!(alpha_all, DataFrame(subject=sub, condition=cond, time=times, power=avg(alpha)))
        for d in ALPHA_CANDIDATES
            append!(alpha_cand[d], DataFrame(subject=sub, condition=cond, time=times,
                                            power=avg(alpha_by[d])))
        end
    end
end

CSV.write(joinpath(RESULTS, "trials.csv"), trials_all)
CSV.write(joinpath(RESULTS, "waveforms.csv"), waves_all)
CSV.write(joinpath(RESULTS, "heog.csv"), heog_all)
CSV.write(joinpath(RESULTS, "qc.csv"), qc_all)
CSV.write(joinpath(RESULTS, "nsw.csv"), nsw_all)
CSV.write(joinpath(RESULTS, "alpha.csv"), alpha_all)
USE_ICA && CSV.write(joinpath(RESULTS, "ica.csv"), ica_all)
for d in ALPHA_CANDIDATES
    CSV.write(joinpath(RESULTS, "alpha_$(d).csv"), alpha_cand[d])
end
serialize(joinpath(RESULTS, "roi_trials.jls"), roi_trials)   # input of script 05
@info "Done: $(length(subjects)) participants, $(count(trials_all.keep)) trials kept."
