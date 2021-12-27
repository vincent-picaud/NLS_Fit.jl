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
plot(X,Y, seriestype = :scatter, label = "raw data")
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

This package embeds the
[NLS_Solver.jl](https://github.com/vincent-picaud/NLS_Solver.jl)
package to perform nonlinear least squares regressions. To use it you
simply have to add the `NLS_Solver` prefix. The example below explains
this in details.

We want to solve this nonlinear least squares problem:
```math
\min\limits_\theta \frac{1}{2}\|Y-m(X,\theta)\|_2^2
```

The objective function is created as follows:

```@example session
nls = NLS_ForwardDiff_From_Model2Fit(model,X,Y)
```

We now select the solver, here the Levenberg-Marquardt method

```@example session
conf = NLS_Solver.LevenbergMarquardt_Conf()
```

and solve the problem:

```@example session
result = NLS_Solver.solve(nls,θ_init,conf)
```

Note that we have used the `NLS_Solver.` prefix, **without** importing
again the `NLS_Solver.jl` package.

!!! danger "Caveat"
    Do not re-import the `NLS_Solver.jl` package by typing `using NLS_Solver`. 
	As `NLS_Fit` already embeds this package, this may cause versions conflicts.
	The right way to go is to add a `NLS_Solver` prefix when you need
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

