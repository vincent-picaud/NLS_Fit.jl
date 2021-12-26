using Revise
using NLS_Fit
using DelimitedFiles

rootDir = joinpath(dirname(pathof(NLS_Fit)), "..")
dataDir = joinpath(rootDir,"data")
dataFile = joinpath(dataDir,"simple_gaussian.txt")

XY=readdlm(dataFile)
X=XY[:,1]
Y=XY[:,2]

model = Gaussian_Peak() 

θ_init = Float64[1,10,5]

nls = NLS_ForwardDiff_From_Model2Fit(model,X,Y)
conf = NLS_Solver.Levenberg_Marquardt_Conf()
