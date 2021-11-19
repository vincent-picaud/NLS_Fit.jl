 ```@meta
CurrentModule = NLS_Fit
```

```@setup session
using NLS_Fit
using DelimitedFiles

using Plots
ENV["GKSwstype"]=100
gr()

rootDir  = joinpath(dirname(pathof(NLS_Fit)), "..")
dataDir = joinpath(rootDir,"data")
```

# Simple fit

**Plot data :**

```@example session
XY=readdlm(joinpath(dataDir,"simple_gaussian.txt")) # hide
X = XY[:,1] # hide
Y = XY[:,2] # hide
plot(X,Y, seriestype = :scatter, label = "raw data", title = "Simple 1D Plot")
```

**Prepare model and initial θ :**

```@example session
model = Gaussian_Peak()
θ_init = Float64[1,10,5]
Y_init = eval_y(model,X,θ_init)
plot!(X,Y_init, label = "model θ_init")
```

**Wrap and call a NLS_Solver :**

```@example session
using NLS_Solver

nls = NLS_ForwardDiff_From_Model2Fit(model,X,Y)
conf = Levenberg_Marquardt_Conf()
result = solve(nls,θ_init,conf)
```

**Use result :**

```@example session
converged(result)
```

```@example session
θ_fit = solution(result)
Y_fit = eval_y(model,X,θ_fit)
plot!(X,Y_fit, label = "model θ_fit")
```

# Fit with recalibration

This problem is 3 Gaussian at positions 5, 10, 20. However the loaded data `(X,Y)` presents a miss calibrated X. The "true" X is:
```math
X = 1.1\ X_\text{true} + 0.2
```

This example use synthetic data that can simulated as follows:

```julia
using DelimitedFiles,Random
using NLS_Fit

Random.seed!(1234)

model = Gaussian_Peak() + Gaussian_Peak() + Gaussian_Peak()
θ1 = Float64[1,5,1]
θ2 = Float64[1,10,1]
θ3 = Float64[1,20,1]
θ = vcat(θ1,θ2,θ3)

X=Float64[1:0.25:30;]
n=length(X)
Y=eval_y(model,X,θ) + 0.1*(rand(n) .- 0.5)

@. X = 1.1*X + 0.2

writedlm("simple_recalibration.txt",hcat(X,Y))
```

**Plot data :**

```@example session
XY=readdlm(joinpath(dataDir,"simple_recalibration.txt")) # hide
X = XY[:,1] # hide
Y = XY[:,2] # hide
plot(X,Y, seriestype = :scatter, label = "raw data", title = "Recalibration")
```


**Prepare model and initial θ :**

We first define the uncalibrated model

```@example session
model = Gaussian_Peak() + Gaussian_Peak() + Gaussian_Peak()

θ1 = Float64[1,5,1]
θ2 = Float64[2,10,1]
θ3 = Float64[1,20,2]
θ_init_model = vcat(θ1,θ2,θ3)
```

Then this model is complete with a parameterized transformation. Here we use a [`Map_Affine_Monotonic`](@ref) with initial parameters `θc = Float64[1,0]`

```@example session
recal_map = Map_Affine_Monotonic(X[1],X[end])
recal_model = Recalibration(model,recal_map)
	
θ_map = Float64[1,0]
θ_init_recal_model = vcat(θ_init_model, θ_map)
	
Y_init = eval_y(recal_model,X,θ_init_recal_model)
	
plot!(X,Y_init, label = "model θ_init")	
```

**Wrap and call a NLS_Solver :**

We must constrain positions, we use a bound constrained solver.

```@example session
ε = eps(Float64)
lower_bound = Float64[0,5,ε,0,10,ε,0,20,ε,0.5,0.0] 
upper_bound = Float64[+Inf,5,2.5,+Inf,10,2.5,+Inf,20,2.5,1.5,2.0]
bc = BoundConstraints(lower_bound,upper_bound)
```

```@example session
nls = NLS_ForwardDiff_From_Model2Fit(recal_model,X,Y)
conf = Levenberg_Marquardt_BC_Conf()
result = solve(nls,θ_init_recal_model,bc,conf)
```

**Use result :**

```@example session
converged(result)
```

```@example session
θ_fit_recal_model = solution(result)
Y_fit_recal_model = eval_y(recal_model,X,θ_fit_recal_model)
plot!(X,Y_fit_recal_model, label = "model θ_fit")
```

**Recalibrate X :**

To recalibrate X we can apply the fitted transformation:

```@example session
X_recal = eval_x(recal_model,X,θ_fit_recal_model)
```

TODO: to fix

To plot fitted model using this recalibrated X, one must use `model`
(and not `recal_model`). Do not forget to pop the to last calibration
parameters from `θ`:

```@example session
Y_recal = eval_y(model,X,@view θ_fit_recal_model[1:NLS_Fit.parameter_size(model)])
plot(X_recal,Y, label = "recalibrated data")
plot!(X_recal,Y_recal, label = "fitted + recalibrated model")
```

# Position dependant parameters

I this example there are 3 Gaussian with a position dependant shape
factor σ.  For this example an affine law is assumed.


**Plot data :**

```@example session
XY=readdlm(joinpath(dataDir,"varying_σ.txt")) # hide
X = XY[:,1] # hide
Y = XY[:,2] # hide
plot(X,Y, seriestype = :scatter, label = "raw data", title = "Gaussians with a position dependant shape
factor")
```

**Prepare model and initial θ :**

We first create a model with 3 Gaussian peaks

```@example session
model = Gaussian_Peak() + Gaussian_Peak() + Gaussian_Peak()

θ1 = Float64[1,5,1]
θ2 = Float64[1,10,1]
θ3 = Float64[1,20,1]

θ_model = vcat(θ1,θ2,θ3)
```

In this model the three shape factors σ1, σ2, σ3 are all equal to
one. We want to create a model where σ follows an affine law:

```math
σ(X) = L_A(X) σ_A θ_A + L_B(X) σ_B θ_B 
```
where ``L_A, L_B`` is the Lagrange basis.

```@example session
# define a map:  X_A => σ_A,  X_B => σ_B
map_pos2sigma = Map_Affine(1.0  => 1.0, 30.0 => 5.0)
# initial parameter value
θ_map = ones(NLS_Fit.parameter_size(map_pos2sigma))
```

We now create a vector with Gaussian positions μ1, μ2, μ3 that will
the positions to map to find σ(X). We also need the indices of the
three σ parameters (that is their indices in θ_model vector).

```@example session
indices = [3,6,9];
pos     = [5.0, 10.0, 20.0];
```

We now have all the required information to create the model with a
varying shape factor:

```@example session
model_with_σ_law = Model2Fit_Mapped_Parameters(model,map_pos2sigma,indices,pos);
```

The new parameter vector `θ_model_with_σ_law` is build from the initial `θ_model`. We must first remove the parameters associated to σ1, σ2, σ3 (as there are going to be replaced by σ(X)) and then add the σ(X) map own parameters:

```@example session
θ_model_with_σ_law = vcat(deleteat!(copy(θ_model),indices),θ_map)
```

There is a helper function [`get_model_θ`](@ref) allowing to retrieve easily the parameters **after** transformation:

```@example session
get_model_θ(model_with_σ_law,θ_model_with_σ_law)
```

We see that our first approximation overestimate greatly σ3. This is
even more obvious when plotting the model:

```@example session
Y_model_with_σ_law = eval_y(model_with_σ_law,X,θ_model_with_σ_law)
plot!(X,Y_model_with_σ_law, label = "model σ law, init")
```

To perform a nonlinear least squares fitting the procedure is as usual:

```@example session
nls = NLS_ForwardDiff_From_Model2Fit(model_with_σ_law,X,Y)
conf = Levenberg_Marquardt_Conf()
result = solve(nls,θ_model_with_σ_law,conf)
```

```@example session
converged(result)
```

```@example session
solution(result)
```

```@example session
θ_fitted = solution(result)
Y_fitted = eval_y(model_with_σ_law,X,θ_fitted)
plot!(X,Y_fitted, label = "model θ_fit")
```

As before, we can get back the fitted parameters of the individual 3 Gaussian peaks:

```@example session
get_model_θ(model_with_σ_law,solution(result))
```

For comparison the true solution is:

```julia
9-element Vector{Float64}:
  1.0
  5.0
  1.0
  1.0
 10.0
  1.5
  1.0
 20.0
  2.0
```
