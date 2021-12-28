# Generate 3 Gaussian peak at position : 5, 10, 20 with difference
# shape factors σ 1, 1.5, 2
#
using DelimitedFiles,Random
using NLS_Fit

Random.seed!(1234)

model = Gaussian_Peak() + Gaussian_Peak() + Gaussian_Peak()
θ1 = Float64[1,5,1]
θ2 = Float64[1,10,1.5]
θ3 = Float64[1,20,2]
θ = vcat(θ1,θ2,θ3)

X=Float64[1:0.25:30;]
n=length(X)
Y=eval_y(model,X,θ) + 0.1*(rand(n) .- 0.5)

rootDir = joinpath(dirname(pathof(NLS_Fit)), "..")
dataDir = joinpath(rootDir,"data")
dataFile = joinpath(dataDir,"constrained_parameters.txt")
writedlm(dataFile,hcat(X,Y))

