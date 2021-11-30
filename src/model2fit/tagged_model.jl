# Tag model and embed data
#
export Model2Fit_TaggedModel
export get_data 

@doc raw"""
```julia
Model2Fit_TaggedModel(model,data)
```

Tag model and embed data

# Extra method

- [`get_data`](@ref) 
"""
struct Model2Fit_TaggedModel{MODEL<:Abstract_Model2Fit,DATA} <: Abstract_Model2Fit
    _model::MODEL
    _data::DATA
end

# Extra methods  ================
#
@doc raw"""
```julia
get_data(m::Model2Fit_TaggedModel{MODEL,DATA)::DATA
```

Return embedded data
"""
get_data(m::Model2Fit_TaggedModel) = m._data
get_model(m::Model2Fit_TaggedModel) = m._model

# Visit  ================
#
visit_submodel_size(model::Model2Fit_TaggedModel) = 1
visit_get_submodel(model::Model2Fit_TaggedModel,submodel_idx::Int) = get_model(model)
visit_get_Y(model::Model2Fit_TaggedModel,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,θ::AbstractVector) = Y
visit_get_X(model::Model2Fit_TaggedModel,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,θ::AbstractVector) = X
visit_get_θ(model::Model2Fit_TaggedModel,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,θ::AbstractVector) = θ


# Interface  ================
#
parameter_size(m::Model2Fit_TaggedModel) = parameter_size(visit_get_submodel(m,1))

function accumulate_y!(m::Model2Fit_TaggedModel,Y::AbstractVector,X::AbstractVector,θ::AbstractVector)
    @assert length(θ) == parameter_size(m)
    @assert length(X) == length(Y)

    accumulate_y!(get_model(model),Y,X,θ)
end

