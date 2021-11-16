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

```@example session
XY=readdlm(joinpath(dataDir,"simple_gaussian.txt")) 
X = XY[:,1]
Y = XY[:,2]
plot(X,Y)
```
