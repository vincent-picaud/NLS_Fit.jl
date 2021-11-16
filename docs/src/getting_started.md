```@setup session
using NLS_Fit
using DelimitedFiles

using Plots
ENV["GKSwstype"]=100
gr()

# ENV["GKSwstype"]=100
# using GR

rootDir  = joinpath(dirname(pathof(NLS_Fit)), "..")
dataDir = joinpath(rootDir,"data")
```

# Simple 1D Plot

```@repl session
Y = rand(10)
plot(Y)
```
