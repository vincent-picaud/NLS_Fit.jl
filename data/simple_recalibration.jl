# Generate 3 Gaussian peak at position : 5, 10, 20 but transform true
# X into 1.1 * X + 0.2 before saving
#
using DelimitedFiles,Random
using NLS_Fit

Random.seed!(1234)

n = 30
model = Gaussian_Peak() + Gaussian_Peak() + Gaussian_Peak()
θ1 = Float64[1,5,1]
θ2 = Float64[1,10,1]
θ3 = Float64[1,20,1]
θ = vcat(θ1,θ2,θ3)

X=Float64[1:n;]
Y=eval_y(model,X,θ) + 0.1*(rand(n) .- 0.5)

@. X = 1.1*X + 0.2

writedlm("simple_recalibration.txt",hcat(X,Y))

