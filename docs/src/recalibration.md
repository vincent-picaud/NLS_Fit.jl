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
θ2 = Float64[1.5,10,2]
θ3 = Float64[0.75,20,3]
θ = vcat(θ1,θ2,θ3)

X=Float64[1:0.25:30;]
n=length(X)
Y=eval_y(model,X,θ) + 0.1*(rand(n) .- 0.5)

@. X = 1.1*X + 0.2 # inverse map is: 0.91 * X - 0.18

rootDir = joinpath(dirname(pathof(NLS_Fit)), "..")
dataDir = joinpath(rootDir,"data")
dataFile = joinpath(dataDir,"recalibration.txt")
writedlm(dataFile,hcat(X,Y))


```

If we plot this, we get:

```@example session
XY=readdlm(joinpath(dataDir,"recalibration.txt")) # hide
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

θ1 = Float64[1.0,5,1]
θ2 = Float64[1.5,10,2]
θ3 = Float64[0.75,20,3]
θ_uncalibrated_model = vcat(θ1,θ2,θ3)
```

This model is the sum of three Gaussian peaks, its parameter vector θ is
```math
\theta = [ h_1, \mu_1, \sigma_1, h_2, \mu_2, \sigma_2, h_3, \mu_3, \sigma_3]
```

Even if the centers ``\mu_1, \mu_2, \mu_3`` are corrects, the peaks are miss positioned due to the bad calibration:

```@example session
Y_uncalibrated_model = eval_y(model,X,θ_uncalibrated_model)

plot(X,Y, seriestype = :scatter, label = "raw data", title = "Recalibration")
plot!(X,Y_uncalibrated_model, label = "uncalibrated model")
```

## The recalibrable model

The idea to create a recalibrable model ``\hat{m}`` is as follows. We use an initial model ``m(\theta,\X)`` to create a new one where the spatial variable ``\hat{X}`` is transformed ``X=f_{\hat{\theta}}(\hat{X})``:
```math
\hat{m}([\theta,\hat{\theta}],\hat{X})=m(\theta,X=f_{\hat{\theta}}(\hat{X}))
```
With this approach, the initial model ``m`` sees the calibrated ``X``. The calibration map, ``f_{\hat{\theta}}`` is parameterized by ``\hat{\theta}``. The create calibrable model is defined by a new parameter vector of the form ``[\theta,\hat{\theta}]``.

For the current example we will use an affine function for ``f_{\hat{\theta}}``. This function is defined as follows:

```math
f_{\hat{\theta}}(\hat{X}) = L_A(\hat{X})X_A\hat{\theta}_A + L_B(\hat{X})X_B\hat{\theta}_B
```
where ``L_A, L_B`` are the Lagrange basis, 
```math
L_A(\hat{X})=\frac{\hat{X}_B-\hat{X}}{\hat{X}_B-\hat{X}_A},\ \ L_B(\hat{X})=\frac{\hat{X}-\hat{X}_A}{\hat{X}_B-\hat{X}_A}
```
We remark that for ``\hat{\theta}= [\hat{\theta}_A,\hat{\theta}_B]=[1,1]`` the transform ``f_{\hat{\theta}}`` is the transform that maps  ``\hat{X}_A`` to ``X_A`` and ``\hat{X}_B`` to ``X_B``. This did not happen by chance.

!!! note 
    Whenever possible we try to use parametrizations such that ``\hat{\theta}=[1,\dots,1]`` has a 
	special role, like identity transform. This ``\hat{\theta}=[1,\dots,1]`` can be used as parameter 
	initial value. This is important from a numerical point of view, where it is better to have an 
	unknown vector ``\theta`` where all components have the same order of magnitude.
   
On the Julia side, the affine transform ``f_{\hat{\theta}}`` is declared as follows:
```julia
Map_Affine(hat_X_A => X_A, hat_X_B => X_B)
```
whereas the shortened version
```julia
Map_Affine(X_A, X_B)
```
means
```julia
Map_Affine(X_A => X_A, X_B => X_B)
```
In that case ``\hat{\theta}=[1,1]`` is the identity transform.

This is the syntax we use to define our ``f_{\hat{\theta}}``:

```@example session
recalibration_map = Map_Affine(X[1],X[end])
```

Now the recalibrable model ``\hat{m}`` is defined as follows:

```@example session
recalibration_model = Model2Fit_Recalibration(model,recalibration_map)
```

Now our unknown vector ``[\theta,\hat{\theta}]`` is:

```@example session
θ_map = Float64[1,1]
θ_init_recalibration_model = vcat(θ_uncalibrated_model, θ_map)
```
which is, in a more readable form, equivalent to:

```math
\theta_{\text{recalibrable}} = [ h_1, \mu_1, \sigma_1, h_2, \mu_2, \sigma_2, h_3, \mu_3, \sigma_3, \hat{\theta}_A, \hat{\theta}_B]
```

## The bound constrained solver













