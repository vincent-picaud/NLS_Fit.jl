using NLS_Fit
using Test, BenchmarkTools

@testset "NLS_Fit.jl" begin

    include("misc/Misc.jl")    

    include("map/Map.jl")    

    include("peak/Peak.jl")    

    include("model2fit_sum.jl")    

    include("solver_wrappers/Solver_Wrappers.jl")    
    
end
