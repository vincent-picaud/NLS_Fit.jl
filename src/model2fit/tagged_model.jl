# Tag model and embed data
#
export Model2Fit_TaggedModel
export get_data, get_tagged_model # TODO rename get_data -> get_tagged_data
export get_tagged_data_type, get_tagged_model_type

@doc raw"""
```julia
Model2Fit_TaggedModel(model,data)
```

Tag model and embed data

# Extra method

- [`get_data`](@ref) 
- [`get_tagged_model`](@ref) 
- [`get_tagged_data_type`](@ref) 
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
get_tagged_model(m::Model2Fit_TaggedModel) = m._model

@doc raw"""
```julia
get_tagged_data_type(m::Abstract_Model2Fit)::DataType
```

Return embedded data type, or `Nothing` if m is not a
[`Model2Fit_TaggedModel`](@ref).
"""
get_tagged_data_type(m::Abstract_Model2Fit) = Nothing
get_tagged_data_type(m::Model2Fit_TaggedModel{MODEL,DATA}) where {MODEL<:Abstract_Model2Fit,DATA} = DATA


@doc raw"""
```julia
get_tagged_model_type(m::Abstract_Model2Fit)::DataType
```

Return wrapped model type, or `Nothing` if m is not a
[`Model2Fit_TaggedModel`](@ref).
"""
get_tagged_model_type(m::Abstract_Model2Fit) = Nothing
get_tagged_model_type(m::Model2Fit_TaggedModel{MODEL,DATA}) where {MODEL<:Abstract_Model2Fit,DATA} = MODEL


# Visit  ================
#
visit_submodel_size(model::Model2Fit_TaggedModel) = 1
visit_get_submodel(model::Model2Fit_TaggedModel,submodel_idx::Int) = get_tagged_model(model)
visit_get_Y(model::Model2Fit_TaggedModel,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,θ::AbstractVector) = Y
visit_get_X(model::Model2Fit_TaggedModel,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,θ::AbstractVector) = X
visit_get_θ(model::Model2Fit_TaggedModel,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,θ::AbstractVector) = θ


# Interface  ================
#
parameter_size(m::Model2Fit_TaggedModel) = parameter_size(visit_get_submodel(m,1))

function accumulate_y!(m::Model2Fit_TaggedModel,Y::AbstractVector,X::AbstractVector,θ::AbstractVector)
    @assert length(θ) == parameter_size(m)
    @assert length(X) == length(Y)

    accumulate_y!(get_tagged_model(m),Y,X,θ)
end

