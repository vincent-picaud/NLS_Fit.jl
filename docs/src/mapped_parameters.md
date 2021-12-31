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

# Mapped parameters

You can reproduce the following computation using `sandbox/mapped_parameters.jl`.

## Data and model

In this example there are 3 Gaussian peaks with different shape factors ``\sigma_1, \sigma_2, \sigma_3``. The goal here is to constrain these shape factors to have an affine dependence with respect to ``X``.

Let us plot the initial data and model.

```@example session
model = Gaussian_Peak() + Gaussian_Peak() + Gaussian_Peak()

θ1 = Float64[1,5,1]
θ2 = Float64[1,10,1]
θ3 = Float64[1,20,1]

θ_model = vcat(θ1,θ2,θ3)
```

This parameter vector is:
```math
[h_1, \mu_1, \sigma_1, h_2, \mu_2, \sigma_2, h_3, \mu_3, \sigma_3]
```

```@example session
XY=readdlm(joinpath(dataDir,"mapped_parameters.txt")) # hide
X = XY[:,1] # hide
Y = XY[:,2] # hide
plot(X,Y, seriestype = :scatter, label = "raw data")

Y_model = eval_y(model,X,θ_model)
plot!(X,Y_model, label = "initial model")
```
For the model the three shape factors are ``\sigma_i=1``. The data has been generated with ``(\sigma_1, \sigma_2, \sigma_3=(1,1.5,2)``.

It often happens that we want to constraint some parameters to follow a given law. Here we want to create a model where σ follows an affine law:

```math
σ(X) = L_A(X) σ_A θ_A + L_B(X) σ_B θ_B 
```
where ``L_A, L_B`` are the Lagrange basis.

A detailed description about this affine function parametrization is
given in [A model for calibration](@ref A_model_for_calibration).

```@example session
# define a map:  X_A => σ_A,  X_B => σ_B
map_pos2sigma = Map_Affine(1.0  => 1.0, 30.0 => 5.0)
# initial parameter value
θ_map = ones(NLS_Fit.parameter_size(map_pos2sigma))
nothing # hide
```

We now create a vector of components ``[\mu_1, \mu_2, \mu_3]``. This
components are the 3 Gaussian peak positions. These positions are used
find ``\sigma_i=\sigma(μ_i)`` using the `map_pos2sigma`. To write
these shape factors at the right emplacements in the `θ_model` vector
we also need the indices of the ``\sigma_1, \sigma_2, \sigma_3``
parameters.

```@example session
σ_indices = [3,6,9]
ref_pos   = [5.0, 10.0, 20.0]
nothing # hide
```

We now have all the required information to create the model where the
3 shape factors follow the affine law ``\sigma(X)``:

```@example session
model_with_σ_law = Model2Fit_Mapped_Parameters(model,map_pos2sigma,σ_indices,ref_pos)
nothing # hide
```

The new parameter vector `θ_model_with_σ_law` is build from the
initial parameter vector `θ_model`. We must first remove the
parameters associated to ``\sigma_1, \sigma_2, \sigma_3`` (as there
are going to be replaced by ``\sigma(\mu_1), \sigma(\mu_2),
\sigma(\mu_3)``). We then have to add the `map_pos2sigma` extra
parameters. All these modifications can be performed by:

```@example session
θ_model_with_σ_law = vcat(deleteat!(copy(θ_model),σ_indices),θ_map)
```

This parameter vector is of the form:
```math
[h_1, \mu_1, h_2, \mu_2, h_3, \mu_3, \theta_A, \theta_B]
```

We can easily retrieve the parameters vector `θ_model` with values induces by the `map_pos2sigma` transformation:

```@example session
get_model_θ(model_with_σ_law,θ_model_with_σ_law)
```

We see that the model with its initial `θ_model_with_σ_law` parameter vector overestimate the shape factor ``\sigma_i``. This is also visible if we plot the model ``Y`` value:

```@example session
Y_model_with_σ_law = eval_y(model_with_σ_law,X,θ_model_with_σ_law)
plot!(X,Y_model_with_σ_law, label = "initial model with constrained σ")
```

## Fit
