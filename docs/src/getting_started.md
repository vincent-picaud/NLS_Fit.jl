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

```@example session
model = Gaussian_Peak() + Gaussian_Peak() + Gaussian_Peak()
recal_model = Recalibration_Affine(model,X[1],X[end])

θ1 = Float64[1,5,1]
θ2 = Float64[1,10,1]
θ3 = Float64[1,20,1]
θ_init_model = vcat(θ1,θ2,θ3)

# calibration parameters
θc = Float64[1,1]
θ_init_recal_model = vcat(θ_init_model, θc)

Y_init = eval_y(recal_model,X,θ_init_recal_model)

plot!(X,Y_init, label = "model θ_init")
```

**Wrap and call a NLS_Solver :**

We must constrain positions, we use a bound constrained solver.

```@example session
using NLS_Solver
ε = eps(Float64)
lower_bound = Float64[0,5,ε,0,10,ε,0,20,ε,0.5,0.5]
upper_bound = Float64[+Inf,5,2.5,+Inf,10,2.5,+Inf,20,2.5,1.5,1.5]
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

To plot fitted model using this recalibrated X, one must use `model`
(and not `recal_model`). Do not forget to pop the to last calibration
parameters from `θ`:

```@example session
Y_recal = eval_y(model,X,@view θ_fit_recal_model[1:NLS_Fit.parameter_size(model)])
plot(X_recal,Y, label = "recalibrated data")
plot!(X_recal,Y_recal, label = "fitted + recalibrated model")
```
