export Model2Fit_Mapped_Parameters
export get_model, get_model_θ

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

parameter_size(mp::Model2Fit_Mapped_Parameters) = parameter_size(mp._model) - length(mp._indices) + parameter_size(mp._map)

@doc raw"""
```julia
get_model(mp::Model2Fit_Mapped_Parameters) -> Absatrct_Model2Fit
```

Get back the wrapped model
"""
get_model(mp::Model2Fit_Mapped_Parameters) = mp._model

@doc raw"""
```julia
get_model_θ(mp::Model2Fit_Mapped_Parameters,θ::AbstractVector) -> θ::AbstractVector
```

Retrieve the parameter vector θ associated to the wrapped model [`get_model`](@ref).
"""
function get_model_θ(mp::Model2Fit_Mapped_Parameters,θ::AbstractVector)
    @assert length(θ) == parameter_size(mp)
    
    θ_model_reduced = @view θ[1:(parameter_size(mp._model) - length(mp._indices))]
    θ_model_reduced_size = length(θ_model_reduced)
    
    θ_map = @view θ[(θ_model_reduced_size+1):end]
    mapped_elements = eval_x(mp._map,mp._elements2map, θ_map)
    
    insert_some_elements(θ_model_reduced, mp._indices,mapped_elements)
end

function eval_y!(mp::Model2Fit_Mapped_Parameters,Y::AbstractVector,X::AbstractVector,θ::AbstractVector)
    θ_model = get_model_θ(mp,θ)
    model = get_model(mp)

    eval_y!(model,Y,X,θ_model)
end