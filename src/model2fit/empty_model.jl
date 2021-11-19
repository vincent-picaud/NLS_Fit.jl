export Model2Fit_Empty

"""

An empty model without parameters. This is useful to initialize
summation:

```julia
model = Model2Fit_Empty()

for i in 1:n
    model = model + GaussianPeak()
end
```
"""
struct Model2Fit_Empty <: Abstract_Model2Fit
end

parameter_size(m::Model2Fit_Empty) = 0
eval_y!(m::Model2Fit_Empty,Y::AbstractVector,X::AbstractVector,θ::AbstractVector) = Y
