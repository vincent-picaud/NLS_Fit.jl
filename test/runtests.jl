using NLS_Fit
using Test, BenchmarkTools

@testset "NLS_Fit.jl" begin

    include("peak/Peak.jl")    

    include("model2fit_sum.jl")    
    
end
