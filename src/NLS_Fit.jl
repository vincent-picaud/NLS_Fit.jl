module NLS_Fit

include("misc/Misc.jl")

include("map/Map.jl")

include("abstract_model2fit.jl")

include("model2fit_empty.jl")
include("model2fit_sum.jl")

include("peak/Peak.jl")

include("recalibration/Recalibration.jl")

include("solver_wrappers/Solver_Wrappers.jl")

end
