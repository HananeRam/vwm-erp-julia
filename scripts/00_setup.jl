# Run once from the project folder:  julia --project=. scripts/00_setup.jl
# Creates Project.toml and Manifest.toml (commit both, so others get the
# same package versions).
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.add(["EDF", "DSP", "CSV", "DataFrames", "MixedModels", "Effects", "Plots",
         "HypothesisTests", "Distributions", "MultivariateStats", "MAT"])
Pkg.instantiate()
println("Setup complete.")
