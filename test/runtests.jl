using SingleFactorANOVA
using Test

@testset "SingleFactorANOVA.jl" begin
    X = [[1, 2, 5, 9], [2, 6, 4, 2, 3, 8], [15, 6, 26]]
    result = anova(X)
    @test typeof(result) == AnovaResult
    @test result.SSB == 303.44230769230774
    @test result.SSE == 268.25
    @test result.DFB == 2
    @test result.DFE == 10
    @test result.MSB == 151.72115384615387
    @test result.MSE == 26.825
    @test result.F == 5.655961000788588
    @test result.p == 0.022745050729447377
    tk_result = tukey_kramer(X, result)
    @test typeof(tk_result, TukeyKramerResult)
    @test getindex(tk_result.mean_differences,1) == 0.08333333333333304
    @test getindex(tk_result.q_crit,1) == 12.458712776723205
    @test getindex(tk_result.mean_differences,1) == false
end
