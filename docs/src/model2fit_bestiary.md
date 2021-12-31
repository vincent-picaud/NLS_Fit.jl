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

# `Model2Fit` bestiary

We try to catalog all the available model here.

## Peak models

### Gaussian peak

See: [`Gaussian_Peak`](@ref)

## Model or Parameter transformations

### Recalibration

See [`Model2Fit_Recalibration`](@ref) .

### Shared parameters

See [`Model2Fit_Shared_Parameters`](@ref) 

### Mapped parameters

See [`Model2Fit_Mapped_Parameters`](@ref) 

### Constant parameters

TODO

### Stacked parameters

Useful when you have several ROIs.

See [`Model2Fit_Stacked`](@ref) 
