using SingleFactorANOVA
using Test

@testset "SingleFactorANOVA.jl" begin
    result = anova([[1, 2, 5, 9], [2, 6, 4, 2, 3, 8], [15, 6, 26]])
    @test typeof(result) == AnovaResult
    @test result.SSB == 303.44230769230774
    @test result.SSE == 268.25
    @test result.DFB == 2
    @test result.DFE == 10
    @test result.MSB == 151.72115384615387
    @test result.MSE == 26.825
    @test result.F == 5.655961000788588
    @test result.p == 0.022745050729447377
end
