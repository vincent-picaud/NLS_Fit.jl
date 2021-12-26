import NLS_Solver

using ForwardDiff

export NLS_ForwardDiff_From_Model2Fit

# See details in NLS_ForwardDiff_From_Model2Fit doc
#
# In short: by export NLS_Solver, we avoid potential version
# compatibility problems. The only inconvenience is that user have to
# use the "NLS_Solver." prefix.
#
export NLS_Solver

# Important: ################
#
# To reduce number of memory allocations it is important to fully
# parameterize the NLS_ForwardDiff_From_Model2Fit structure
#
# NLS_ForwardDiff_From_Model2Fit{MODEL2FIT_TYPE <: Abstract_Model2Fit,
#                                X_ELEMENT_TYPE,
#                                Y_ELEMENT_TYPE,
#                                X_TYPE <: AbstractVector{X_ELEMENT_TYPE},
#                                Y_TYPE <: AbstractVector{Y_ELEMENT_TYPE}} <: NLS_Solver.AbstractNLS
#
# - Without parametrization ( NLS_ForwardDiff_From_Model2Fit <: NLS_Solver.AbstractNLS ):
#
# julia> @btime eval_r($nls,$θ)
#   2.724 μs (44 allocations: 928 bytes)
#
# - With parametrization:
#
# julia> @btime eval_r($nls,$θ)
#   842.590 ns (1 allocation: 160 bytes)
#
@doc raw"""

A wrapper that allows to use the
[NLS_Solver.jl](https://github.com/vincent-picaud/NLS_Solver.jl)
package to solve the nonlinear least squares problem associated to a
model plus its (X,Y) data.

You must construct an instance of nls problem as follows:
```julia
NLS_ForwardDiff_From_Model2Fit(fit_model::MODEL2FIT_TYPE,
                               X::X_TYPE,
                               Y::Y_TYPE)
```

You can then use the `NLS_Solver` package as usual, but with prefixed
`NLS_Solver.`

```@example
nls = NLS_ForwardDiff_From_Model2Fit(model,X,Y)
conf = NLS_Solver.LevenbergMarquardt_Conf()
result = NLS_Solver.solve(nls,θ_init,conf)
```

**CAVEAT:**

Please note that you should **not** explicitly
```@example
using NLS_Solver
```
as this package is reexported from `NLS_Fit`. 
Doing so avoids potential version compatibility problems.

"""
struct NLS_ForwardDiff_From_Model2Fit{MODEL2FIT_TYPE <: Abstract_Model2Fit,
                                      X_ELEMENT_TYPE,
                                      Y_ELEMENT_TYPE,
                                      X_TYPE <: AbstractVector{X_ELEMENT_TYPE},
                                      Y_TYPE <: AbstractVector{Y_ELEMENT_TYPE}} <: NLS_Solver.AbstractNLS
    _fit_model::MODEL2FIT_TYPE
    _X::X_TYPE
    _Y::Y_TYPE
 
    function NLS_ForwardDiff_From_Model2Fit(fit_model::MODEL2FIT_TYPE,
                                            X::X_TYPE,Y::Y_TYPE) where{MODEL2FIT_TYPE <: Abstract_Model2Fit,
                                                                       X_ELEMENT_TYPE,
                                                                       Y_ELEMENT_TYPE,
                                                                       X_TYPE <: AbstractVector{X_ELEMENT_TYPE},
                                                                       Y_TYPE <: AbstractVector{Y_ELEMENT_TYPE}} 
        @assert length(X) == length(Y)

        new{MODEL2FIT_TYPE,X_ELEMENT_TYPE,Y_ELEMENT_TYPE,X_TYPE,Y_TYPE}(fit_model,X,Y)
    end
end



NLS_Solver.parameter_size(nls::NLS_ForwardDiff_From_Model2Fit) = parameter_size(nls._fit_model)
NLS_Solver.residue_size(nls::NLS_ForwardDiff_From_Model2Fit)  = length(nls._Y)

function NLS_Solver.eval_r(nls::NLS_ForwardDiff_From_Model2Fit,θ::AbstractVector)
    Y_model = eval_y(nls._fit_model,nls._X,θ)
    @. Y_model = nls._Y - Y_model
    Y_model
end


function NLS_Solver.eval_r_J(nls::NLS_ForwardDiff_From_Model2Fit, θ::AbstractVector)
    r = NLS_Solver.eval_r(nls,θ)
    J = ForwardDiff.jacobian(θ->NLS_Solver.eval_r(nls,θ), θ)
    r,J
end

