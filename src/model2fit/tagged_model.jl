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

# Visit  ================
#
visit_submodel_size(model::Model2Fit_TaggedModel) = 1

_visit_get_submodel(model::Model2Fit_TaggedModel) = model._model

function visit_get_submodel(model::Model2Fit_TaggedModel,submodel_idx::Int)
    @assert 1 ≤ submodel_idx ≤ visit_submodel_size(model)

    _visit_get_submodel(model)
end

_visit_get_X(model::Model2Fit_TaggedModel,X::AbstractVector) = X

function visit_get_X(model::Model2Fit_TaggedModel,submodel_idx::Int,X::AbstractVector,θ::AbstractVector)
    @assert 1 ≤ submodel_idx ≤ visit_submodel_size(model)
    @assert length(θ) == parameter_size(model)
    
    _visit_get_X(model,X)
end

_visit_get_θ(model::Model2Fit_TaggedModel,θ::AbstractVector) = θ
function visit_get_θ(model::Model2Fit_TaggedModel,submodel_idx::Int,X::AbstractVector,θ::AbstractVector)
    @assert 1 ≤ submodel_idx ≤ visit_submodel_size(model)
    @assert length(θ) == parameter_size(model)

    _visit_get_θ(model,θ)
end

# Interface  ================
#
parameter_size(m::Model2Fit_TaggedModel) = parameter_size(_visit_get_submodel(m))

function accumulate_y!(m::Model2Fit_TaggedModel,Y::AbstractVector,X::AbstractVector,θ::AbstractVector)
    @assert length(θ) == parameter_size(m)
    @assert length(X) == length(Y)

    accumulate_y!(_visit_get_submodel(m),Y,_visit_get_X(m,X),_visit_get_θ(m,θ))

    Y
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
