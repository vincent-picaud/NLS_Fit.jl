# Create a model with const parameters ****************
#
export Model2Fit_Const_Parameters
export get_model, get_model_θ

@doc raw"""

Set some parameters to constant values

# Example

Let's assume that we have an initial model `model` with `θ` as parameter vector. 
If we want θ1, θ3, θ5 to be constant, proceed as follows

```julia
indices = [1,3,5]
values = Float64[1,2,3]
        
model_with_const_params = Model2Fit_Const_Parameters(model, indices, values)
```

Do not forget to udpate parameter vector θ, this can be done as
follows:

```julia
deleteat!(θ, indices)
```

# Also see

"""
struct Model2Fit_Const_Parameters{MODEL <: Abstract_Model2Fit,
                                  INDICES <: AbstractVector{Int},
                                  VALUES <: AbstractVector}  <: Abstract_Model2Fit
    _model::MODEL
    _indices::INDICES
    _values::VALUES

    function Model2Fit_Const_Parameters(model::MODEL,
                                        indices::INDICES,
                                        values::VALUES) where {MODEL <: Abstract_Model2Fit,
                                                               INDICES <: AbstractVector{Int},
                                                               VALUES <: AbstractVector}
        @assert (length(indices) == 0) || ( (1 ≤ minimum(indices)) && (maximum(indices) ≤ parameter_size(model)) )
        @assert allunique(indices)
        @assert length(indices) == length(values)
        
        new{MODEL,INDICES,VALUES}(model, indices, values)

    end
    
end

# Convenience methods ================
#
@doc raw"""
```julia
get_model(model::Model2Fit_Const_Parameters) -> Absatrct_Model2Fit
```

Get back the wrapped model
"""
get_model(model::Model2Fit_Const_Parameters) = model._model

@doc raw"""
```julia
get_model_θ(mp::Model2Fit_Const_Parameters,θ::AbstractVector) -> θ::AbstractVector
```

Retrieve the parameter vector θ associated to the wrapped model [`get_model`](@ref).
"""
function get_model_θ(model::Model2Fit_Const_Parameters,θ::AbstractVector)
    @assert length(θ) == parameter_size(model)
    
    # Create the new θ by inserting constant values
    #
    insert_some_elements(θ, model._indices, model._values)
end

# Visit  ================
#
# note: no need to check
#   @assert 1 ≤ submodel_idx ≤ visit_submodel_size(model)
#   @assert length(X) == length(Y)
# as this is done higher in the visit() method
#
visit_submodel_size(model::Model2Fit_Const_Parameters) = 1
visit_get_submodel(model::Model2Fit_Const_Parameters,submodel_idx::Int) = model._model
visit_get_X(model::Model2Fit_Const_Parameters,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,θ::AbstractVector) = X
visit_get_Y(model::Model2Fit_Const_Parameters,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,θ::AbstractVector) = Y
visit_get_θ(model::Model2Fit_Const_Parameters,submodel_idx::Int,Y::AbstractVector,X::AbstractVector,θ::AbstractVector) = get_model_θ(model,θ)


# Interface ================
#
parameter_size(mp::Model2Fit_Const_Parameters) = parameter_size(mp._model) - length(mp._indices) 

function accumulate_y!(model::Model2Fit_Const_Parameters,Y::AbstractVector,X::AbstractVector,θ::AbstractVector)
    @assert length(θ) == parameter_size(model)
    
    inner_model = get_model(model)
    θ_inner_model = get_model_θ(model,θ)

    accumulate_y!(inner_model,Y,X,θ_inner_model)
end

