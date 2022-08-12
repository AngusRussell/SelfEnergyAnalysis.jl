module SelfEnergyAnalysis

# Write your package code here.
using CSV
using DataFrames
using Glob
using Hilbert
using Dierckx
using LsqFit
using Statistics
using LaTeXStrings
using PyPlot
#pygui(true)

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
