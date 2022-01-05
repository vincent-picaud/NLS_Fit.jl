#
# Recursively with model
#
# Important: this is the tool one must use to extract back information
# after fitting. Instead of trying to track individual parameters,
# simply visit the models and extract required information.
#
# Note that you can also tag the models that interest you using a
# Model2Fit_TaggedModel in which you can embed some extra data.
#
export visit

# ================================================================
# visit model specializations
# ================================================================
#
@doc raw"""
```julia
visit(mp::Abstract_Model2Fit,X::AbstractVector,θp::AbstractVector,action::Function) -> nothing
```

Recursively visit models using a depth-first-order.

For each visited model, performian `action` which is a function of
type

```julia
visit_default_action(m::Abstract_Model2Fit,x::AbstractVector,θ::AbstractVector)  -> Bool
```

If the function returns `false` the depth-first-search is stopped.

# Implementation details

The `visit` functionality requires these methods to be defined for each visited model:
- `visit_submodel_size(model)`
- `visit_get_submodel(model,submodel_idx)`
- `visit_get_Y(model,submodel_idx,Y,X,θ)`
- `visit_get_X(model,submodel_idx,Y,X,θ)`
- `visit_get_θ(model,submodel_idx,Y,X,θ)`

There is no need to export these methods.

These functions **never** modifies `Y` components, in some cases like
[`Model2Fit_Stacked`](@ref), `visit_get_Y` returns a view.
"""
function
visit(action::Function,mp::Abstract_Model2Fit,Y::AbstractVector,X::AbstractVector,θ::AbstractVector)
    # as we perform a sanity check there is no need to repeat it
    # in all the visit_... functions
    @assert length(θ) == parameter_size(mp)
    @assert length(X) == length(Y)
    
    continue_exploration = action(mp,Y,X,θ)
    
    if continue_exploration
        n_submodel = visit_submodel_size(mp)
        for submodel_idx in 1:n_submodel
           
            visit(action,
                  visit_get_submodel(mp,submodel_idx),
                  visit_get_Y(mp,submodel_idx,Y,X,θ),
                  visit_get_X(mp,submodel_idx,Y,X,θ),
                  visit_get_θ(mp,submodel_idx,Y,X,θ))
        end
    end
end

# Convience methods ================
#
function visit_debug(mp::Abstract_Model2Fit,Y::AbstractVector,X::AbstractVector,θ::AbstractVector)
    visit(mp,Y,X,θ) do m::Abstract_Model2Fit,y::AbstractVector,x::AbstractVector,θ::AbstractVector
        @assert parameter_size(m)==size(θ,1)
        println("model type :",typeof(m))
        println("parameters :",θ)
        true
    end
end
