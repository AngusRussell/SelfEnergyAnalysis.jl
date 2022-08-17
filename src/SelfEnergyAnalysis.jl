module SelfEnergyAnalysis

# Write your package code here.
using CSV   # https://github.com/JuliaData/CSV.jl
using DataFrames    # https://github.com/JuliaData/DataFrames.jl 
using Glob  # https://github.com/JuliaPackageMirrors/Glob.jl
using Hilbert   # https://github.com/fpreiswerk/Hilbert.jl
using Dierckx   # https://github.com/kbarbary/Dierckx.jl
using LsqFit    # https://github.com/JuliaNLSolvers/LsqFit.jl
using Statistics    # https://github.com/JuliaStats/Statistics.jl
using LaTeXStrings  # https://github.com/stevengj/LaTeXStrings.jl
using PyPlot    # https://github.com/JuliaPy/PyPlot.jl
using Markdown

include("plotting_presets.jl")
include("process_raw_data.jl")
include("useful_manipulations.jl") 
include("MDC.jl") 
include("EDC.jl") 
include("self_energy.jl") 
include("lower_polariton_dispsion.jl") 
include("spectral_function.jl") 
include("bare-band_algorithm.jl") 

export fig, load_all, crop_data, attenuation_correction_scalar, remove_background, correct_BK30, correct_BK40, MDC, Plot_MDC_Cut_Max, Plot_MDC_Cut_Fit, energy_band_fit, energy_band_max, ImagExtraction, RealExtraction, Σ_E, SpectralFunction


end
