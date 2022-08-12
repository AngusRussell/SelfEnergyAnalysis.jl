using SelfEnergyAnalysis
using Test



@testset "SelfEnergyAnalysis.jl" begin
    # Write your tests here.
    a = [1, NaN, 2, NaN, 3]
    SelfEnergyAnalysis.remove_NaN!(a)
    @test any(isnan.(a)) == false
    b = [ NaN, NaN, NaN, NaN]
    SelfEnergyAnalysis.remove_NaN!(b)
    @test any(isnan.(b)) == false
end
