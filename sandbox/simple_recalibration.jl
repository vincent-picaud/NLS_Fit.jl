using BenchmarkTools
using Revise
using NLS_Fit
using NLS_Solver
using DelimitedFiles

XY=readdlm("../data/simple_recalibration.txt")
X=XY[:,1]
Y=XY[:,2]

model = Gaussian_Peak() + Gaussian_Peak() + Gaussian_Peak()
model = Recalibration_Affine(model,X[1],X[end])

θ1 = Float64[1,5,1]
θ2 = Float64[1,10,1]
θ3 = Float64[1,20,1]
# calibration parameters
θc = Float64[1,1]

θ = vcat(θ1,θ2,θ3,θc)

nls = NLS_ForwardDiff_From_Model2Fit(model,X,Y)
conf = Levenberg_Marquardt_Conf()
result = solve(nls,θ,conf)

eval_x(model,X,solution(result))
eval_x(model,X,θ)




@benchmark eval_x($model,$X,$θ)

