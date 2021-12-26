var documenterSearchIndex = {"docs":
[{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"CurrentModule = NLS_Fit","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"using NLS_Fit\nusing DelimitedFiles\n\nusing Plots\nENV[\"GKSwstype\"]=100\ngr()\n\nrootDir  = joinpath(dirname(pathof(NLS_Fit)), \"..\")\ndataDir = joinpath(rootDir,\"data\")","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"We present here some examples of increasing complexity.","category":"page"},{"location":"getting_started/#Simple-fit","page":"Simple fit","title":"Simple fit","text":"","category":"section"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"The first example is a simple Gaussian peak fit. We use the data/simple_gaussian.txt data file.","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"XY=readdlm(joinpath(dataDir,\"simple_gaussian.txt\")) # hide\nX = XY[:,1] # hide\nY = XY[:,2] # hide\nplot(X,Y, seriestype = :scatter, label = \"raw data\", title = \"Simple 1D Plot\")","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"The model here is really simple, a single Gaussian peak. One must also provide a parameter vector θ. For this model, θ=[h,μ,σ] where h is peak height, μ its center and σ its shape factor (see Gaussian_Peak for further details). The function eval_y() eval model Y values given X values and its parameter vector θ.","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"model = Gaussian_Peak()\nθ_init = Float64[1,10,5]\nY_init = eval_y(model,X,θ_init)\nplot!(X,Y_init, label = \"initial model\")","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"Wrap and call a NLS_Solver :","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"using NLS_Solver\n\nnls = NLS_ForwardDiff_From_Model2Fit(model,X,Y)\nconf = Levenberg_Marquardt_Conf()\nresult = solve(nls,θ_init,conf)","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"Use result :","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"converged(result)","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"θ_fit = solution(result)\nY_fit = eval_y(model,X,θ_fit)\nplot!(X,Y_fit, label = \"fitted model\")","category":"page"},{"location":"getting_started/#Fit-with-recalibration","page":"Simple fit","title":"Fit with recalibration","text":"","category":"section"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"This problem is 3 Gaussian peaks at positions 5, 10, 20. However the loaded data (X,Y) presents a miscalibrated X. This miscalibrated X is computed from The \"true\" X as follows:","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"X = 11 X_texttrue + 02","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"The complete process to generate the synthetic data is:","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"using DelimitedFiles,Random\nusing NLS_Fit\n\nRandom.seed!(1234)\n\nmodel = Gaussian_Peak() + Gaussian_Peak() + Gaussian_Peak()\nθ1 = Float64[1,5,1]\nθ2 = Float64[1,10,1]\nθ3 = Float64[1,20,1]\nθ = vcat(θ1,θ2,θ3)\n\nX=Float64[1:0.25:30;]\nn=length(X)\nY=eval_y(model,X,θ) + 0.1*(rand(n) .- 0.5)\n\n@. X = 1.1*X + 0.2\n\nwritedlm(\"simple_recalibration.txt\",hcat(X,Y))","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"If we plot this, we get:","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"XY=readdlm(joinpath(dataDir,\"simple_recalibration.txt\")) # hide\nX = XY[:,1] # hide\nY = XY[:,2] # hide\nplot(X,Y, seriestype = :scatter, label = \"raw data\", title = \"Recalibration\")","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"We can also plot miscalibrated data and the \"true\" one that can be retrieved by inverting the X(X_texttrue) relation:","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"X_true = (X .- 0.2)/1.1\nplot(X_true, seriestype = :scatter, label = \"true X\", title = \"X-axis miscalibration\")\nplot!(X, seriestype = :scatter, label = \"miscalibrated X\")","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"The first step is to define an uncalibrated model.","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"model = Gaussian_Peak() + Gaussian_Peak() + Gaussian_Peak()\n\nθ1 = Float64[1,5,1]\nθ2 = Float64[2,10,1]\nθ3 = Float64[1,20,2]\nθ_uncalibrated_model = vcat(θ1,θ2,θ3)","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"Then this model is complete with a parameterized transformation. Here we use a Map_Affine_Monotonic with initial parameters θ_map = Float64[θ_A,θ_B]. The advantage of such transformation is that we can impose an increasing map with simple bound constraints θ_B > 1. More precisely, by Map_Affine_Monotonic(X[1],X[end]) we define an identity recalibration map when θ_map = [1.0, 0.0]","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"recalibration_map = Map_Affine_Monotonic(X[1],X[end])\nrecalibration_model = Model2Fit_Recalibration(model,recalibration_map)\n\t\nθ_map = Float64[1,1]\nθ_init_recalibration_model = vcat(θ_uncalibrated_model, θ_map)\n\t\nY_init = eval_y(recalibration_model,X,θ_init_recalibration_model)\n\t\nplot(X,Y_init, label = \"initial model\")\t","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"Wrap and call a NLS_Solver :","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"We must constrain positions, we use a bound constrained solver.","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"ε = eps(Float64)\nlower_bound = Float64[0,5,ε,0,10,ε,0,20,ε,0.5,0.0] \nupper_bound = Float64[+Inf,5,2.5,+Inf,10,2.5,+Inf,20,2.5,1.5,2.0]\nbc = BoundConstraints(lower_bound,upper_bound)","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"nls = NLS_ForwardDiff_From_Model2Fit(recalibration_model,X,Y)\nconf = Levenberg_Marquardt_BC_Conf()\nresult = solve(nls,θ_init_recalibration_model,bc,conf)","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"Use result :","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"converged(result)","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"θ_fit_recalibration_model = solution(result)\nY_fit_recalibration_model = eval_y(recalibration_model,X,θ_fit_recalibration_model)\nplot!(X,Y_fit_recalibration_model, label = \"fitted model\")","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"Recalibrate X :","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"To recalibrate X we can apply the fitted transformation:","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"X_recal = eval_calibrated_x(recalibration_model,X,θ_fit_recalibration_model)","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"TODO: to fix","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"To plot fitted model using this recalibrated X, one must use model (and not recalibration_model). Do not forget to pop the to last calibration parameters from θ:","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"Y_recal = eval_y(model,X,@view θ_fit_recalibration_model[1:NLS_Fit.parameter_size(model)])\nplot(X_recal,Y, label = \"recalibrated data\")\nplot!(X_recal,Y_recal, label = \"fitted + recalibrated model\")","category":"page"},{"location":"getting_started/#Position-dependant-parameters","page":"Simple fit","title":"Position dependant parameters","text":"","category":"section"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"I this example there are 3 Gaussian with a position dependant shape factor σ.  For this example an affine law is assumed.","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"Plot data :","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"XY=readdlm(joinpath(dataDir,\"varying_σ.txt\")) # hide\nX = XY[:,1] # hide\nY = XY[:,2] # hide\nplot(X,Y, seriestype = :scatter, label = \"raw data\", title = \"Gaussians with a position dependant shape\nfactor\")","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"Prepare model and initial θ :","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"We first create a model with 3 Gaussian peaks","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"model = Gaussian_Peak() + Gaussian_Peak() + Gaussian_Peak()\n\nθ1 = Float64[1,5,1]\nθ2 = Float64[1,10,1]\nθ3 = Float64[1,20,1]\n\nθ_model = vcat(θ1,θ2,θ3)","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"In this model the three shape factors σ1, σ2, σ3 are all equal to one. We want to create a model where σ follows an affine law:","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"σ(X) = L_A(X) σ_A θ_A + L_B(X) σ_B θ_B ","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"where L_A L_B is the Lagrange basis.","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"# define a map:  X_A => σ_A,  X_B => σ_B\nmap_pos2sigma = Map_Affine(1.0  => 1.0, 30.0 => 5.0)\n# initial parameter value\nθ_map = ones(NLS_Fit.parameter_size(map_pos2sigma))","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"We now create a vector with Gaussian positions μ1, μ2, μ3 that will the positions to map to find σ(X). We also need the indices of the three σ parameters (that is their indices in θ_model vector).","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"indices = [3,6,9];\npos     = [5.0, 10.0, 20.0];","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"We now have all the required information to create the model with a varying shape factor:","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"model_with_σ_law = Model2Fit_Mapped_Parameters(model,map_pos2sigma,indices,pos);","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"The new parameter vector θ_model_with_σ_law is build from the initial θ_model. We must first remove the parameters associated to σ1, σ2, σ3 (as there are going to be replaced by σ(X)) and then add the σ(X) map own parameters:","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"θ_model_with_σ_law = vcat(deleteat!(copy(θ_model),indices),θ_map)","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"There is a helper function get_model_θ allowing to retrieve easily the parameters after transformation:","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"get_model_θ(model_with_σ_law,θ_model_with_σ_law)","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"We see that our first approximation overestimate greatly σ3. This is even more obvious when plotting the model:","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"Y_model_with_σ_law = eval_y(model_with_σ_law,X,θ_model_with_σ_law)\nplot!(X,Y_model_with_σ_law, label = \"initial model\")","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"To perform a nonlinear least squares fitting the procedure is as usual:","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"nls = NLS_ForwardDiff_From_Model2Fit(model_with_σ_law,X,Y)\nconf = Levenberg_Marquardt_Conf()\nresult = solve(nls,θ_model_with_σ_law,conf)","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"converged(result)","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"solution(result)","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"θ_fitted = solution(result)\nY_fitted = eval_y(model_with_σ_law,X,θ_fitted)\nplot!(X,Y_fitted, label = \"fitted model\")","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"As before, we can get back the fitted parameters of the individual 3 Gaussian peaks:","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"get_model_θ(model_with_σ_law,solution(result))","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"For comparison the true solution is:","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"9-element Vector{Float64}:\n  1.0\n  5.0\n  1.0\n  1.0\n 10.0\n  1.5\n  1.0\n 20.0\n  2.0","category":"page"},{"location":"getting_started/#Position-dependant-parameters-2","page":"Simple fit","title":"Position dependant parameters","text":"","category":"section"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"I this example there are 3 Gaussian with a position dependant shape factor σ.  For this example an affine law is assumed.","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"Plot data :","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"XY=readdlm(joinpath(dataDir,\"varying_σ.txt\")) # hide\nX = XY[:,1] # hide\nY = XY[:,2] # hide\nplot(X,Y, seriestype = :scatter, label = \"raw data\", title = \"Gaussians with a position dependant shape\nfactor\")","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"Prepare model and initial θ :","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"We first create a model with 3 Gaussian peaks","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"model = Gaussian_Peak() + Gaussian_Peak() + Gaussian_Peak()\n\nθ1 = Float64[1,5,1]\nθ2 = Float64[1,10,1]\nθ3 = Float64[1,20,1]\n\nθ_model = vcat(θ1,θ2,θ3)","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"In this model the three shape factors σ1, σ2, σ3 are all equal to one. We want to create a model where σ follows an affine law:","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"σ(X) = L_A(X) σ_A θ_A + L_B(X) σ_B θ_B ","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"where L_A L_B is the Lagrange basis.","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"# define a map:  X_A => σ_A,  X_B => σ_B\nmap_pos2sigma = Map_Affine(1.0  => 1.0, 30.0 => 5.0)\n# initial parameter value\nθ_map = ones(NLS_Fit.parameter_size(map_pos2sigma))","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"We now create a vector with Gaussian positions μ1, μ2, μ3 that will the positions to map to find σ(X). We also need the indices of the three σ parameters (that is their indices in θ_model vector).","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"σ_indices = [3,6,9];\nref_pos     = [5.0, 10.0, 20.0];","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"We now have all the required information to create the model with a varying shape factor:","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"model_with_σ_law = Model2Fit_Mapped_Parameters(model,map_pos2sigma,σ_indices,ref_pos);","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"The new parameter vector θ_model_with_σ_law is build from the initial θ_model. We must first remove the parameters associated to σ1, σ2, σ3 (as there are going to be replaced by σ(X)) and then add the σ(X) map own parameters:","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"θ_model_with_σ_law = vcat(deleteat!(copy(θ_model),σ_indices),θ_map)","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"There is a helper function get_model_θ allowing to retrieve easily the parameters after transformation:","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"get_model_θ(model_with_σ_law,θ_model_with_σ_law)","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"We see that our first approximation overestimate greatly σ3. This is even more obvious when plotting the model:","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"Y_model_with_σ_law = eval_y(model_with_σ_law,X,θ_model_with_σ_law)\nplot!(X,Y_model_with_σ_law, label = \"initial model\")","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"To perform a nonlinear least squares fitting the procedure is as usual:","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"nls = NLS_ForwardDiff_From_Model2Fit(model_with_σ_law,X,Y)\nconf = Levenberg_Marquardt_Conf()\nresult = solve(nls,θ_model_with_σ_law,conf)","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"converged(result)","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"solution(result)","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"θ_fitted = solution(result)\nY_fitted = eval_y(model_with_σ_law,X,θ_fitted)\nplot!(X,Y_fitted, label = \"fitted model\")","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"As before, we can get back the fitted parameters of the individual 3 Gaussian peaks:","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"get_model_θ(model_with_σ_law,solution(result))","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"For comparison the true solution is:","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"9-element Vector{Float64}:\n  1.0\n  5.0\n  1.0\n  1.0\n 10.0\n  1.5\n  1.0\n 20.0\n  2.0","category":"page"},{"location":"getting_started/#Position-dependant-parameters-and-recalibration","page":"Simple fit","title":"Position dependant parameters and recalibration","text":"","category":"section"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"In this example we add a recalibration. We proceed as before.","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"XY=readdlm(joinpath(dataDir,\"varying_σ_and_recalibration.txt\")) # hide\nX = XY[:,1] # hide\nY = XY[:,2] # hide\n\nrecalibration_map = Map_Affine_Monotonic(X[1],X[end])\nmodel_with_σ_law_and_recal = Model2Fit_Recalibration(model_with_σ_law,recalibration_map)\n\nθ_map = Float64[1,1]\nθ_model_with_σ_law_and_recal = vcat(θ_model_with_σ_law, θ_map)","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"The new uncalibrated data and the initial model are:","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"\nY_init = eval_y(model_with_σ_law_and_recal,X,θ_model_with_σ_law_and_recal)\nplot(X,Y, seriestype = :scatter, label = \"raw data\")\nplot!(X,Y_init, label = \"initial model\")\n","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"We solve the problem as before","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"nls = NLS_ForwardDiff_From_Model2Fit(model_with_σ_law_and_recal,X,Y)\nconf = Levenberg_Marquardt_Conf()\nresult = solve(nls,θ_model_with_σ_law_and_recal,conf)","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"and plot the solution","category":"page"},{"location":"getting_started/","page":"Simple fit","title":"Simple fit","text":"θ_fitted = solution(result)\nY_fitted = eval_y(model_with_σ_law_and_recal,X,θ_fitted)\n\nplot(X,Y, seriestype = :scatter, label = \"raw data\")\nplot!(X,Y_fitted, label = \"fitted model\")\n","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = NLS_Fit","category":"page"},{"location":"#NLS_Fit","page":"Home","title":"NLS_Fit","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for NLS_Fit.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [NLS_Fit]","category":"page"},{"location":"#NLS_Fit.Abstract_Map","page":"Home","title":"NLS_Fit.Abstract_Map","text":"Map base type. \n\nThis type is an abstraction of a fX  Y=f(X).\n\nSuch transformations can be used in various contexts:\n\ncalibration task: we have a reference hatX to a calibrated X\nparameter transformation: we map a parameter vector hatθ to a new one θ.\n\nInterface\n\nparameter_size\neval_map \n\n\n\n\n\n","category":"type"},{"location":"#NLS_Fit.Abstract_Model2Fit","page":"Home","title":"NLS_Fit.Abstract_Model2Fit","text":"Abstract type, base of all model to fit.\n\nInterface\n\nparameter_size \naccumulate_y! \n\n\n\n\n\n","category":"type"},{"location":"#NLS_Fit.Gaussian_Peak","page":"Home","title":"NLS_Fit.Gaussian_Peak","text":"Gaussian peak\n\np(x θ = hμσ) = h e^-frac12 left(fracx-μσ right)^2\n\n\n\n\n\n","category":"type"},{"location":"#NLS_Fit.Map_Affine","page":"Home","title":"NLS_Fit.Map_Affine","text":"Map_Affine(X_hat_A => X_A, X_hat_B => X_B)\n\nMap_Affine(X_hat_A, X_hat_B)\n\nThe second constructor is a shortcut for\n\nMap_Affine(X_hat_A => X_hat_A, X_hat_B => X_hat_B)\n\nDefines an affine transformation:\n\nX(hatX) = L_A(hatX) X_A θ_A +  L_B(hatX) X_B θ_B\n\nwhere L_A and L_B are the Lagrange basis.\n\nIn peculiar, for (θ_Aθ_B) = (11) we have an affine map such that X(hatX_A) = X_A and X(hatX_B) = X_B.\n\n\n\n\n\n","category":"type"},{"location":"#NLS_Fit.Map_Affine_Monotonic","page":"Home","title":"NLS_Fit.Map_Affine_Monotonic","text":"Map_Affine_Monotonic(X_hat_A => X_A, X_hat_B => X_B)\n\nMap_Affine_Monotonic(X_hat_A, X_hat_B)\n\nThe second constructor is a shortcut for\n\nMap_Affine_Monotonic(X_hat_A => X_hat_A, X_hat_B => X_hat_B)\n\nSame as Map_Affine but uses another parametrization that allows to insure monotonic map using simple bound constraints.\n\nThe Map_Affine map is defined as follows:\n\nX(hatX) = L_A(hatX) X_A θ_A +  L_B(hatX) X_B θ_B\n\nWith this parametrization you have to add the X_B θ_B  X_A θ_A to insure that the map is increasing.\n\nOn the other hand, with the Map_Affine_Monotonic parametrization:\n\nX(hatX) = L_A(hatX) X_A θ_A +  L_B(hatX) ( (X_B-X_A) θ_B + X_A θ_A )\n\nwhere \n\nfracdXdhatX = fracX_B-X_AhatX_B-hatX_A θ_B\n\na simply bound constraint θ_B  0 is sufficient to preserve the increasing or decreasing character. This preserved slope is given by the sign of fracX_B-X_AhatX_B-hatX_A which is constant once the structure has been initialized.\n\nIn peculiar, for (θ_Aθ_B) = (11) we have an affine map such that X(hatX_A) = X_A and X(hatX_B) = X_B. Whereas for (θ_Aθ_B) = (10) we have a constant one.\n\n\n\n\n\n","category":"type"},{"location":"#NLS_Fit.Map_Translate","page":"Home","title":"NLS_Fit.Map_Translate","text":"Map_Translate(;scale = 1)\n\nDefines a translate transformation:\n\nX(hatX) = hatX + scale (θ - 1)\n\nIn peculiar, for (θ) = (1) the transformation is identity.  This is due to the presence of the θ- 1 factor. This factor is introduced to stay consistent with other transformations where we try to have identity transform for parameter vector of ones.\n\n\n\n\n\n","category":"type"},{"location":"#NLS_Fit.Model2Fit_Empty","page":"Home","title":"NLS_Fit.Model2Fit_Empty","text":"An empty model without parameters. This is useful to initialize summation:\n\nmodel = Model2Fit_Empty()\n\nfor i in 1:n\n    model = model + GaussianPeak()\nend\n\n\n\n\n\n","category":"type"},{"location":"#NLS_Fit.Model2Fit_Mapped_Parameters","page":"Home","title":"NLS_Fit.Model2Fit_Mapped_Parameters","text":"Create a new model, where some parameters are computed using a Abstract_Map.\n\nTODO: not clear how to use it: certainly need a refactoring...\n\nAlso see\n\nModel2Fit_Shared_Parameters \n\n\n\n\n\n","category":"type"},{"location":"#NLS_Fit.Model2Fit_Recalibration","page":"Home","title":"NLS_Fit.Model2Fit_Recalibration","text":"Create a calibrable model\n\nConstructors\n\ncalibrable_model = Model2Fit_Recalibration(model_to_calibrate, calibration_map)\n\nwhere model_to_calibrate is an Abstract_Model2Fit and calibration_map is an Abstract_Map\n\nExplanation\n\nThe call\n\neval_y(calibrable_model,X_hat,[θ, θ_map])\n\nis equivalent to \n\neval_y(model_to_calibrate,eval_calibrated_x(X_hat, θ_map), [θ])\n\n\n\n\n\n","category":"type"},{"location":"#NLS_Fit.Model2Fit_Shared_Parameters","page":"Home","title":"NLS_Fit.Model2Fit_Shared_Parameters","text":"Share a set of parameters\n\nExample\n\nLet's assume that we have an initial model model with θ as parameter vector.  If we want θ1, θ3, θ5 to share a same value θshared then:\n\nindices_to_share = [1,3,5]\n        \nmodel_with_shared_params = Model2Fit_Shared_Parameters(model, indices_to_share)\n\nDo not forget to udpate parameter vector θ, this can be done as follows:\n\ndeleteat!(θ, indices_to_share)\npush!(θ, θshared)\n\nAlso see\n\nModel2Fit_Mapped_Parameters \n\n\n\n\n\n","category":"type"},{"location":"#NLS_Fit.Model2Fit_Stacked","page":"Home","title":"NLS_Fit.Model2Fit_Stacked","text":"Define a composite model, where each submodel is responsible of a given ROI\n\n| (1)        | (2)       | ...  |––––––|–––––-| ... | models[1]  | model[2]  | ...\n\nDomain (1) range is defined as 1:ROIsizes[1] Domain (2) range is defined as (ROIsizes[1]+1):ROI_sizes[2] etc...\n\n\n\n\n\n","category":"type"},{"location":"#NLS_Fit.Model2Fit_TaggedModel","page":"Home","title":"NLS_Fit.Model2Fit_TaggedModel","text":"Model2Fit_TaggedModel(model,data)\n\nTag model and embed data\n\nExtra method\n\nget_tagged_data \nget_tagged_model \nget_tagged_data_type \n\n\n\n\n\n","category":"type"},{"location":"#NLS_Fit.accumulate_y!-Tuple{NLS_Fit.Abstract_Model2Fit, AbstractVector, AbstractVector, AbstractVector}","page":"Home","title":"NLS_Fit.accumulate_y!","text":"accumulate_y!(::Abstract_Model2Fit,Y::AbstractVector,X::AbstractVector,θ::AbstractVector)::AbstractVector\n\nAccumulate model contribution into vector Y.\n\nAlso see : \n\nAbstract_Model2Fit\neval_y \n\n\n\n\n\n","category":"method"},{"location":"#NLS_Fit.eval_calibrated_x-Tuple{Model2Fit_Recalibration, AbstractVector, AbstractVector}","page":"Home","title":"NLS_Fit.eval_calibrated_x","text":"eval_calibrated_x(m::Model2Fit_Recalibration,X_hat::AbstractVector,θ::AbstractVector) -> AbstractVector\n\nCompute the calibrated X from the reference domain hatX for the given transformation parameters θ.\n\n\n\n\n\n","category":"method"},{"location":"#NLS_Fit.eval_map-Tuple{NLS_Fit.Abstract_Map, AbstractVector, AbstractVector}","page":"Home","title":"NLS_Fit.eval_map","text":"eval_map(m::Abstract_Map,X_hat::AbstractVector,θ::AbstractVector) -> X::AbstractVector\n\nCompute X=X(hatX).\n\n\n\n\n\n","category":"method"},{"location":"#NLS_Fit.eval_y-Tuple{NLS_Fit.Abstract_Model2Fit, AbstractVector, AbstractVector}","page":"Home","title":"NLS_Fit.eval_y","text":"eval_y(m::Abstract_Model2Fit,X::AbstractVector,θ::AbstractVector)::AbstractVector\n\nA convenience function that call accumulate_y! using a zero initialized Y vector. This returned vector contains model values.\n\nAlso see : \n\nAbstract_Model2Fit\n\n\n\n\n\n","category":"method"},{"location":"#NLS_Fit.get_model-Tuple{Model2Fit_Mapped_Parameters}","page":"Home","title":"NLS_Fit.get_model","text":"get_model(mp::Model2Fit_Mapped_Parameters) -> Absatrct_Model2Fit\n\nGet back the wrapped model\n\n\n\n\n\n","category":"method"},{"location":"#NLS_Fit.get_model-Tuple{Model2Fit_Shared_Parameters}","page":"Home","title":"NLS_Fit.get_model","text":"get_model(mp::Model2Fit_Shared_Parameters) -> Absatrct_Model2Fit\n\nGet back the wrapped model\n\n\n\n\n\n","category":"method"},{"location":"#NLS_Fit.get_model_θ-Tuple{Model2Fit_Mapped_Parameters, AbstractVector}","page":"Home","title":"NLS_Fit.get_model_θ","text":"get_model_θ(mp::Model2Fit_Mapped_Parameters,θ::AbstractVector) -> θ::AbstractVector\n\nRetrieve the parameter vector θ associated to the wrapped model get_model.\n\n\n\n\n\n","category":"method"},{"location":"#NLS_Fit.get_model_θ-Tuple{Model2Fit_Shared_Parameters, AbstractVector}","page":"Home","title":"NLS_Fit.get_model_θ","text":"get_model_θ(mp::Model2Fit_Shared_Parameters,θ::AbstractVector) -> θ::AbstractVector\n\nRetrieve the parameter vector θ associated to the wrapped model get_model.\n\n\n\n\n\n","category":"method"},{"location":"#NLS_Fit.get_tagged_data-Tuple{Model2Fit_TaggedModel}","page":"Home","title":"NLS_Fit.get_tagged_data","text":"get_tagged_data(m::Model2Fit_TaggedModel{MODEL,DATA)::DATA\n\nReturn embedded data\n\n\n\n\n\n","category":"method"},{"location":"#NLS_Fit.get_tagged_data_type-Tuple{NLS_Fit.Abstract_Model2Fit}","page":"Home","title":"NLS_Fit.get_tagged_data_type","text":"get_tagged_data_type(m::Abstract_Model2Fit)::DataType\n\nReturn embedded data type, or Nothing if m is not a Model2Fit_TaggedModel.\n\n\n\n\n\n","category":"method"},{"location":"#NLS_Fit.get_tagged_model-Tuple{Model2Fit_TaggedModel}","page":"Home","title":"NLS_Fit.get_tagged_model","text":"get_tagged_model(m::Model2Fit_TaggedModel{MODEL,DATA)::DATA\n\nReturn tagged model\n\n\n\n\n\n","category":"method"},{"location":"#NLS_Fit.get_tagged_model_type-Tuple{NLS_Fit.Abstract_Model2Fit}","page":"Home","title":"NLS_Fit.get_tagged_model_type","text":"get_tagged_model_type(m::Abstract_Model2Fit)::DataType\n\nReturn wrapped model type, or Nothing if m is not a Model2Fit_TaggedModel.\n\n\n\n\n\n","category":"method"},{"location":"#NLS_Fit.insert_some_elements","page":"Home","title":"NLS_Fit.insert_some_elements","text":"insert_some_elements(X::AbstractVector{T},\n                     indices::AbstractVector{Int},\n                     elements::AbstractVector) -> AbstractVector{T}\n\nInsert elements in X at positions indices.\n\nNote: indices must be strictly ordered.\n\nExample: this functions does a work inverse to the deleteat!\n\njulia> X = Any[1:10;]\n10-element Vector{Any}:\n  1\n  2\n  3\n  4\n  5\n  6\n  7\n  8\n  9\n 10\n\njulia> indices = Int[2, 3, 5, 6, 9, 10];\n\njulia> elements = Float64[2, 3, 5, 6, 9, 10];\n\njulia> deleteat!(X, indices)\n4-element Vector{Any}:\n 1\n 4\n 7\n 8\n\njulia> insert_some_elements(X,indices,elements)\n10-element Vector{Any}:\n  1\n  2.0\n  3.0\n  4\n  5.0\n  6.0\n  7\n  8\n  9.0\n 10.0\n\n\n\n\n\n\n","category":"function"},{"location":"#NLS_Fit.parameter_size-Tuple{NLS_Fit.Abstract_Map}","page":"Home","title":"NLS_Fit.parameter_size","text":"parameter_size(::Abstract_Map) -> Int\n\nReturn θ parameter vector length\n\n\n\n\n\n","category":"method"},{"location":"#NLS_Fit.parameter_size-Tuple{NLS_Fit.Abstract_Model2Fit}","page":"Home","title":"NLS_Fit.parameter_size","text":"parameter_size(::Abstract_Model2Fit)::Int\n\nReturn length of the expected parameter vector θ\n\nAlso see : \n\nAbstract_Model2Fit\n\n\n\n\n\n","category":"method"},{"location":"#NLS_Fit.solve_linear_parameters-Tuple{NLS_Fit.Abstract_Model2Fit, AbstractVector, AbstractVector, AbstractVector, AbstractVector{Int64}}","page":"Home","title":"NLS_Fit.solve_linear_parameters","text":"Find optimal linear parameters\n\n\n\n\n\n","category":"method"},{"location":"#NLS_Fit.visit-Tuple{Function, NLS_Fit.Abstract_Model2Fit, AbstractVector, AbstractVector, AbstractVector}","page":"Home","title":"NLS_Fit.visit","text":"visit(mp::Abstract_Model2Fit,X::AbstractVector,θp::AbstractVector,action::Function) -> nothing\n\nRecursively visit models using a depth-first-order.\n\nFor each visited model, performian action which is a function of type\n\nvisit_default_action(m::Abstract_Model2Fit,x::AbstractVector,θ::AbstractVector)  -> Bool\n\nIf the function returns false the depth-first-search is stopped.\n\nImplementation details\n\nThe visit functionality requires these methods to be defined for each visited model:\n\nvisitsubmodelsize(model)\nvisitgetsubmodel(model,submodel_idx)\nvisitgetY(model,submodel_idx,Y,X,θ)\nvisitgetX(model,submodel_idx,Y,X,θ)\nvisitgetθ(model,submodel_idx,Y,X,θ)\n\nThere is no need to export these methods.\n\nThese functions never modifies Y components, in some cases like Model2Fit_Stacked, visit_get_Y returns a view.\n\n\n\n\n\n","category":"method"}]
}
