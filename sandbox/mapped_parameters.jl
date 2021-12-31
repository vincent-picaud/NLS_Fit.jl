using Revise
using NLS_Fit
using DelimitedFiles
using LinearAlgebra

rootDir = joinpath(dirname(pathof(NLS_Fit)), "..")
dataDir = joinpath(rootDir,"data")
dataFile = joinpath(dataDir,"recalibration.txt")

# Getting started example ****************
#

# Data ================
#
XY=readdlm(dataFile);
X=XY[:,1];
Y=XY[:,2];
