# Data check: sub-01 is stored twice in ds006979, as EDF (read with Beacon's
# EDF.jl) and as BrainVision (read with this project's reader). This script
# confirms both give the same signals. EDF stores 16-bit integers, so small
# rounding differences are expected.
#
# Run:  julia --project=. scripts/03_check_formats.jl

include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "..", "src", "VWMPipeline.jl"))
using .VWMPipeline
using Statistics, Printf

sub = "sub-01"
base = joinpath(BIDS_ROOT, sub, "eeg", "$(sub)_task-$(TASK)_eeg")

edf, fs_edf = load_eeg(base * ".edf"; channels=LOAD_CHANNELS)
bv,  fs_bv  = load_eeg(base * ".vhdr"; channels=LOAD_CHANNELS)

println("Sampling rate  EDF: $fs_edf Hz   BrainVision: $fs_bv Hz")
for ch in LOAD_CHANNELS
    if !(haskey(edf, ch) && haskey(bv, ch))
        println(rpad(ch, 14), " missing in one format")
        continue
    end
    n = min(length(edf[ch]), length(bv[ch]))
    a, b = edf[ch][1:n], bv[ch][1:n]
    @printf("%-14s r = %.6f   max |diff| = %.3f µV   samples EDF/BV = %d/%d\n",
            ch, cor(a, b), maximum(abs.(a .- b)), length(edf[ch]), length(bv[ch]))
end
