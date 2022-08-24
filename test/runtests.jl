using SelfEnergyAnalysis
using Test
using PyPlot
using DataFrames
using CSV
pygui(true)
plt.close("all")
BK30_raw, λ30 = load_all("test\\test_data\\raw_data", "csv", 400, 1340, 3);
BK30_cropped, λ30 = crop_data(BK30_raw, λ30,1);

Plot_all(BK30_raw)

BK30, λ30 = correct_BK30("test\\test_data\\raw_data")
Ebm, kpoints = energy_band_fit(MirrorLeft(BK30_cropped)[:,:,1], λ30, false)
# Plot_MDC_Cut_Max(BK30[:,:,6], 300)
Plot_all_bands(BK30, λ30, "fit", 12)

Eb, kpoints = energy_band_fit(BK30[:,:,1], λ30, true)
Em, kpoints = energy_band_fit(BK30[:,:,8], λ30, true)
# Imaginary self-energy
ImEx, ReKK, LBk, UBk = ImagExtraction(BK30, 8, λ30, -3.8e6, -0.4e6, 0.4e6, 3.8e6);
# Real self-energy
ReEx, ImKK, kpoints = RealExtraction(BK30, 8, λ30, LBk, UBk);
P8_SelfEnergy_df = DataFrame("ReEX"=>ReEx, "ImEx"=>ImEx, "ReKK"=>ReKK, "ImKK"=>ImKK, "kpoints"=>kpoints)
check_symmetry(BK30[:,:,1], λ30)
# Plot_MDC_Cut_Fit(BK30[:,:,6], 300)
# SpectralFunction(ReEx, ImEx, Eb, Em, kpoints)

function BK30_self_energy(pumppower)
    Eb, kpoints = energy_band_fit(BK30[:,:,1], λ30, true)
    Em, kpoints = energy_band_fit(BK30[:,:,pumppower], λ30, true)
    ImEx, ReKK, LBk, UBk =  ImagExtraction(BK30, pumppower, λ30)
    ReEx, ImKK, kpoints = RealExtraction(BK30, pumppower, λ30, LBk, UBk)
    Σ_E(BK30, λ30, ReEx, ImEx, ReKK, ImKK, kpoints)
end

# @testset "SelfEnergyAnalysis.jl" begin
#     # Write your tests here.
#     a = [1, NaN, 2, NaN, 3]
#     SelfEnergyAnalysis.remove_NaN!(a)
#     @test any(isnan.(a)) == false
#     b = [ NaN, NaN, NaN, NaN]
#     SelfEnergyAnalysis.remove_NaN!(b)
#     @test any(isnan.(b)) == false
# end
