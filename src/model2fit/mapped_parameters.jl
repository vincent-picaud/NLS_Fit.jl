export Model2Fit_Mapped_Parameters
export get_model_θ, get_model

@doc raw"""

Create a new model, where some parameters are computed using a
[`Abstract_Map`](@ref).

TODO: not clear how to use it: certainly need a refactoring...

# Also see
- [`Model2Fit_Shared_Parameters`](@ref) 
"""
struct Model2Fit_Mapped_Parameters{MODEL <: Abstract_Model2Fit,
                                   MAP <: Abstract_Map,
                                   INDICES <: AbstractVector{Int},
                                   EL2MAP <: AbstractVector} <: Abstract_Model2Fit
    _model::MODEL
    _map::MAP
    _indices::INDICES
    _elements2map::EL2MAP

    function Model2Fit_Mapped_Parameters(model::MODEL,
                                         map::MAP,
                                         indices::INDICES,
                                         elements2map::EL2MAP) where {MODEL <: Abstract_Model2Fit,
                                                                      MAP <: Abstract_Map,
                                                                      INDICES <: AbstractVector{Int},
                                                                      EL2MAP <: AbstractVector}

        @assert parameter_size(model) ≥ length(indices)
        @assert length(indices) == length(elements2map)

        new{MODEL,MAP,INDICES,EL2MAP}(model,
                                      map,
                                      indices,
                                      elements2map)
    end
    
end

# Visit  ================
#
visit_submodel_size(model::Model2Fit_Mapped_Parameters) = 1

function _visit_get_submodel(model::Model2Fit_Mapped_Parameters)
    model._model
end
function visit_get_submodel(model::Model2Fit_Mapped_Parameters,submodel_idx::Int)
    @assert 1 ≤ submodel_idx ≤ visit_submodel_size(model)

    _visit_get_submodel(model)
end

function _visit_get_X(model::Model2Fit_Mapped_Parameters,X::AbstractVector)
    X
end 
function visit_get_X(model::Model2Fit_Mapped_Parameters,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,θ::AbstractVector)
    _visit_get_X(model,X)
end

function _visit_get_Y(model::Model2Fit_Mapped_Parameters,Y::AbstractVector)
    Y
end 
function visit_get_Y(model::Model2Fit_Mapped_Parameters,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,θ::AbstractVector)
    _visit_get_Y(model,Y)
end 

function _visit_get_θ(model::Model2Fit_Mapped_Parameters,θ::AbstractVector)
    θ_model_reduced = @view θ[1:(parameter_size(model._model) - length(model._indices))]
    θ_model_reduced_size = length(θ_model_reduced)
    
    θ_map = @view θ[(θ_model_reduced_size+1):end]
    mapped_elements = eval_map(model._map,model._elements2map, θ_map)
    
    insert_some_elements(θ_model_reduced, model._indices,mapped_elements)
end

function visit_get_θ(model::Model2Fit_Mapped_Parameters,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,θ::AbstractVector)
    _visit_get_θ(model,θ)
end 


# Interface ================
#
parameter_size(mp::Model2Fit_Mapped_Parameters) = parameter_size(mp._model) - length(mp._indices) + parameter_size(mp._map)

function accumulate_y!(mp::Model2Fit_Mapped_Parameters,Y::AbstractVector,X::AbstractVector,θ::AbstractVector)
    θ_model = _visit_get_θ(mp,θ)
    model = _visit_get_submodel(mp)

    accumulate_y!(model,Y,X,θ_model)
end

# Convenience methods ================
#
@doc raw"""
```julia
get_model(mp::Model2Fit_Mapped_Parameters) -> Absatrct_Model2Fit
```

Get back the wrapped model
"""
get_model(mp::Model2Fit_Mapped_Parameters) = _visit_get_submodel(mp)

@doc raw"""
```julia
get_model_θ(mp::Model2Fit_Mapped_Parameters,θ::AbstractVector) -> θ::AbstractVector
```

Retrieve the parameter vector θ associated to the wrapped model [`get_model`](@ref).
"""
get_model_θ(mp::Model2Fit_Mapped_Parameters,θ::AbstractVector) = _visit_get_θ(mp,θ)
