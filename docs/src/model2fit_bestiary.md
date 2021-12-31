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

We try to catalog all the available 1D models.

## Peak models

### Gaussian peak

See: [`Gaussian_Peak`](@ref)

## Isotopic motifs

Group several peaks together to form a peak motif.

## Baseline models

TODO

## Model or Parameter transformations

### Recalibration

Useful to recalibrate signal knowing some peaks of well defined
positions.


See [`Model2Fit_Recalibration`](@ref) .

### Shared parameters

See [`Model2Fit_Shared_Parameters`](@ref) 

### Mapped parameters

See [`Model2Fit_Mapped_Parameters`](@ref) 

### Constant parameters

See [`Model2Fit_Const_Parameters`](@ref) 

### Stacked parameters

Useful when you have several ROIs.

See [`Model2Fit_Stacked`](@ref) 
