using NLS_Fit
using Test, BenchmarkTools

@testset "NLS_Fit.jl" begin

    include("misc/Misc.jl")    

    include("map/Map.jl")    

    include("model2fit/Model2Fit.jl")    

end
