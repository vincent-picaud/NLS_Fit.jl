 ```@meta
CurrentModule = NLS_Fit
```

```@setup session
#using NLS_Fit
#using DelimitedFiles

using Plots
ENV["GKSwstype"]=100
gr()

#rootDir  = joinpath(dirname(pathof(NLS_Fit)), "..")
#dataDir = joinpath(rootDir,"data")
```

# Simple 1D Plot

```@example session
Y = rand(10)
plot(rand(10),rand(10))
```
