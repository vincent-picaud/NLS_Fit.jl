@testset "model2fit/" begin

    include("model_sum.jl")    
    include("const_parameters.jl")
    
    include("peak/Peak.jl")    

    include("solver_wrappers/Solver_Wrappers.jl")    
end
