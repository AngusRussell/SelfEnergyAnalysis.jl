using SelfEnergyAnalysis
using Test
using PyPlot
pygui(true)

#BK30_raw, λ30 = load_all("test\\test_data\\raw_data", "csv", 400, 1340, 3);
#BK30_cropped, λ30 = crop_data(BK30_raw, λ30, 5);

#SelfEnergyAnalysis.MDC_protocol_fitting(BK30_cropped[:,:,1])
@testset "SelfEnergyAnalysis.jl" begin
    # Write your tests here.
    a = [1, NaN, 2, NaN, 3]
    SelfEnergyAnalysis.remove_NaN!(a)
    @test any(isnan.(a)) == false
    b = [ NaN, NaN, NaN, NaN]
    SelfEnergyAnalysis.remove_NaN!(b)
    @test any(isnan.(b)) == false
end
