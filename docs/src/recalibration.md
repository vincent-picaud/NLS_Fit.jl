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

# Peak fitting and recalibration

## Problem description and data

In this problem one must fit 3 Gaussian peaks at positions 5, 10, 20. However the
loaded data `(X,Y)` presents a miscalibrated X. This miscalibrated X
is computed from The "true" X as follows:

```math
X = 1.1\ X_\text{true} + 0.2
```
The complete process to generate the synthetic data is:

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

If we plot this, we get:

```@example session
XY=readdlm(joinpath(dataDir,"simple_recalibration.txt")) # hide
X = XY[:,1] # hide
Y = XY[:,2] # hide
plot(X,Y, seriestype = :scatter, label = "raw data", title = "Recalibration")
```

We can also plot miscalibrated ``X`` and the "true" one ``X_\text{true}``. We only now ``X``, but as we generated artificial data, we know that 

```math
X = 1.1\ X_\text{true} + 0.2
```

which can easily be inverted to give 


```math
X_\text{true} = (X - 0.2)/1.1
```

We can now plot both ``X`` and ``X_\text{true}``:

```@example session
X_true = (X .- 0.2)/1.1
plot(X_true, seriestype = :scatter, label = "true X", title = "X-axis miscalibration")
plot!(X, seriestype = :scatter, label = "miscalibrated X")
```

## The model

The first step is to define an uncalibrated model.

```@example session
model = Gaussian_Peak() + Gaussian_Peak() + Gaussian_Peak()

θ1 = Float64[1,5,1]
θ2 = Float64[2,10,1]
θ3 = Float64[1,20,2]
θ_uncalibrated_model = vcat(θ1,θ2,θ3)
```

This model is the sum of three Gaussian peaks, its parameter vector θ is
```math
\theta = [ h_1, \mu_1, \sigma_1, h_2, \mu_2, \sigma_2, h_3, \mu_3, \sigma_3]
```

Even if the centers ``\mu_1, \mu_2, \mu_3`` are corrects, the peaks are miss positioned due to the bad calibration:

```@example session
Y_uncalibrated_model = eval_y(model,X,θ_uncalibrated_model

```

