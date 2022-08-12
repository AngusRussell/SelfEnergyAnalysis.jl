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


include("plotting_presets.jl") # add mirror
include("process_raw_data.jl") # complete
include("useful_manipulations.jl") # complete
include("MDC.jl") # complete
include("EDC.jl") # complete
include("self_energy.jl") # Not done 
include("lower_polariton_dispsion.jl") # Not done
include("spectral_function.jl") # Not done
include("bare-band_algorithm.jl") # Not done

export fig, load_all, crop_data, attenuation_correction_scalar, remove_background, correct_BK30, correct_BK40, MDC, Plot_MDC_Cut, energy_band_fit, energy_band_max, ImagExtraction, RealExtraction, Σ_E


end
