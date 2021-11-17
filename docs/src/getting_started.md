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

# Simple 1D Plot

Plot data:

```@example session
XY=readdlm(joinpath(dataDir,"simple_gaussian.txt")) # hide
X = XY[:,1] # hide
Y = XY[:,2] # hide
plot(X,Y, seriestype = :scatter, label = "raw data", title = "Simple 1D Plot")
```

Prepare model and initial θ :

```@example session
model = Gaussian_Peak()
θ_init = Float64[1,10,1]
Y_init = eval_y(model,X,θ_init)
plot!(X,Y_init, label = "model θ_init")
```

Wrap and call a NLS_Solver :

```@example session
using NLS_Solver

nls = NLS_ForwardDiff_From_Model2Fit(model,X,Y)
conf = Levenberg_Marquardt_Conf()
result = solve(nls,θ_init,conf)
```

Use result : 

```@example session
converged(result)
```

```@example session
θ_fit = solution(result)
Y_fit = eval_y(model,X,θ_fit)
plot!(X,Y_fit, label = "model θ_fit")
```
