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

You can reproduce the following computation using
[sandbox/simple_gaussian.jl](../../sandbox/simple_gaussian.jl).

## The data and model
The first example is a simple Gaussian peak fit. We use the
`data/simple_gaussian.txt` data file. 

```@example session
XY=readdlm(joinpath(dataDir,"simple_gaussian.txt")) # hide
X = XY[:,1] # hide
Y = XY[:,2] # hide
plot(X,Y, seriestype = :scatter, label = "raw data", title = "Simple 1D Plot")
```

The model is a single Gaussian peak. 
```@example session
model = Gaussian_Peak()
```
The model does no embeds its model parameters. So before being able to evaluate model values one must define a parameter vector θ. For [`Gaussian_Peak`](@ref), this vector is `[h,μ,σ]`. It stores the peak height, its center and its shape factor.

```@example session
θ_init = Float64[1,10,5]
```

We now have everything we need to evaluate model values. This is the
role of the [`eval_y`](@ref) function. This function compute model Y
values given its parameter vector θ and the evaluation sites X.

```@example session
Y_init = eval_y(model,X,θ_init)
plot!(X,Y_init, label = "initial model")
```

We see that the parameter vector θ is not a good guess to fit the peak
to the raw data. We will see now how to automatically adjust these
parameters using a Levenberg-Marquardt like method.

## Model fitting

This package uses the
[NLS_Solver](https://github.com/vincent-picaud/NLS_Solver.jl) to
perform nonlinear least squares regressions. Please note that, to avoid
potential version problems, it is better **not** `using NLS_Solver` but
using the embedded version that is present in the `NLS_Fit` package. Doing so, you simply have to add the `NLS_Solver` prefix:

```@example session
nls = NLS_ForwardDiff_From_Model2Fit(model,X,Y)
conf = NLS_Solver.LevenbergMarquardt_Conf()
result = NLS_Solver.solve(nls,θ_init,conf)
```

!!! danger "Caveat"
    Do not use `using NLS_Solver` as `NLS_Fit` already embeds the `NLS_Solver.jl` package. 
	To avoid version conflicts, you must use this embedded package version.
	As shown before, you just have to add a `NLS_Solver` prefix when you need
    some of the `NLS_Solver.jl` package functionalities. 
	
## Fit result

We can check that the solver converged

```@example session
NLS_Solver.converged(result)
```
and plot the fitted model

```@example session
θ_fit = NLS_Solver.solution(result)
Y_fit = eval_y(model,X,θ_fit)
plot!(X,Y_fit, label = "fitted model", linewidth=3)
```

