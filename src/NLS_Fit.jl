module NLS_Fit

include("abstract_model2fit.jl")

include("model2fit_empty.jl")
include("model2fit_sum.jl")

include("peak/Peak.jl")

include("solver_wrappers/Solver_Wrappers.jl")

end
