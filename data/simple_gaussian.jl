# Generate some synthetic data
#
using DelimitedFiles,Random
using NLS_Fit

Random.seed!(1234)

n = 50
model = Gaussian_Peak()
θ = Float64[2,n/2,3]

X=Float64[1:n;]
Y=eval_y(model,X,θ) + 0.1*(rand(n) .- 0.5)

writedlm("simple_gaussian.txt",hcat(X,Y))

