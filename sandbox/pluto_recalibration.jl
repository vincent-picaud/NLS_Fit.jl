### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ f497faa1-3b8e-4891-b0e5-679fd29500d7
using Pkg

# ╔═╡ 91ab0997-63de-4ed0-af8d-2221568bef7c
Pkg.activate("..")

# ╔═╡ 0f496f40-14bc-45f7-a890-dbe531e2500d
using NLS_Fit, Plots, DelimitedFiles, NLS_Solver

# ╔═╡ c703884c-de01-4539-b8a2-9fb133d0cad6
md"# Environment & packages"

# ╔═╡ c798b6b4-dd3d-4c33-ab73-8e95eac97613
md"# Raw data plot"

# ╔═╡ 2d2d35af-2e14-495f-a6c6-dc348a60fc2b
begin
	XY = readdlm("../data/simple_recalibration.txt")
	X = XY[:,1]
	Y = XY[:,2] 
	plot(X,Y, seriestype = :scatter, label = "raw data", title = "Recalibration")
end

# ╔═╡ 4dcbd1a8-b79a-43a1-9a8a-6fc5c373dc73
md"# Create model"

# ╔═╡ 29ccb74f-8a62-478d-a849-5610cb2e62c5
md"## Uncalibrable model"

# ╔═╡ 44a9850b-4ecc-42c6-8270-2ffad6903b09
begin
	model = Gaussian_Peak() + Gaussian_Peak() + Gaussian_Peak()
	
	θ1 = Float64[1,5,1]
	θ2 = Float64[2,10,1]
	θ3 = Float64[1,20,2]
	θ_init_model = vcat(θ1,θ2,θ3)
end

# ╔═╡ 21353220-0e6b-40e4-8930-4889f9bb757a
md"## Calibrable model

Complete the previous model with a parameterized transformation used as calibration map
"

# ╔═╡ 5ceea856-87f0-4874-860d-57cf1f31c5e8
begin
	recal_map = Map_Affine_Monotonic(X[1],X[end])
	recal_model = Model2Fit_Recalibration(model,recal_map)

	θc = Float64[1,0]
	θ_init_recal_model = vcat(θ_init_model, θc)

	Y_init = eval_y(recal_model,X,θ_init_recal_model)

	plot!(X,Y_init, label = "model θ_init")	
end

# ╔═╡ fed24429-8e20-4601-8bf8-5a4b401b2534
 NLS_Fit.eval_calibrated_x(recal_model,X,θ_init_recal_model)

# ╔═╡ f331ea35-4d56-4ced-a76e-de621de25593
md"Define parameter bounds"

# ╔═╡ fcc0e63e-d873-4af8-bfe2-55b41e2b9f1a
begin
	ε = eps(Float64)
	lower_bound = Float64[0,5,ε,0,10,ε,0,20,ε,0.5,0.0]
	upper_bound = Float64[+Inf,5,2.5,+Inf,10,2.5,+Inf,20,2.5,1.5,2.0]
	bc = BoundConstraints(lower_bound,upper_bound)
end

# ╔═╡ e61704f1-5eae-4a35-9d5f-0f734e46ef03
md"# Solve the nonlinear least squares problem"

# ╔═╡ e33fe5dd-db88-4de9-8816-b4c36d6765d1
begin
	nls = NLS_ForwardDiff_From_Model2Fit(recal_model,X,Y)
	conf = Levenberg_Marquardt_BC_Conf()
	result = solve(nls,θ_init_recal_model,bc,conf)
end

# ╔═╡ b2e98c3f-5de5-4c23-bbb7-80d39cf219fd
md"Converged?"

# ╔═╡ fb6f9b74-1687-4846-b070-58581a1b234b
converged(result)

# ╔═╡ 4bb856a1-5696-42c8-8e63-fd7abe0eb276
md"Objective function value"

# ╔═╡ 9fa86072-9483-4b63-8dc9-e94469c1a8e9
objective_value(result)

# ╔═╡ 3367e5bd-3159-4327-90a4-01b7799e7064
solution(result)

# ╔═╡ c30df4ad-de11-446c-84ed-33755a78e5aa
md"# Plot calibrated model"

# ╔═╡ e7284fe8-79b8-4684-922c-f004cb98333e
begin
	θ_fit_recal_model = solution(result)
	Y_fit_recal_model = eval_y(recal_model,X,θ_fit_recal_model)
	plot!(X,Y_fit_recal_model, label = "model θ_fit")
end

# ╔═╡ Cell order:
# ╟─c703884c-de01-4539-b8a2-9fb133d0cad6
# ╠═f497faa1-3b8e-4891-b0e5-679fd29500d7
# ╠═91ab0997-63de-4ed0-af8d-2221568bef7c
# ╠═0f496f40-14bc-45f7-a890-dbe531e2500d
# ╟─c798b6b4-dd3d-4c33-ab73-8e95eac97613
# ╠═2d2d35af-2e14-495f-a6c6-dc348a60fc2b
# ╠═4dcbd1a8-b79a-43a1-9a8a-6fc5c373dc73
# ╠═29ccb74f-8a62-478d-a849-5610cb2e62c5
# ╠═44a9850b-4ecc-42c6-8270-2ffad6903b09
# ╟─21353220-0e6b-40e4-8930-4889f9bb757a
# ╠═5ceea856-87f0-4874-860d-57cf1f31c5e8
# ╠═fed24429-8e20-4601-8bf8-5a4b401b2534
# ╠═f331ea35-4d56-4ced-a76e-de621de25593
# ╠═fcc0e63e-d873-4af8-bfe2-55b41e2b9f1a
# ╠═e61704f1-5eae-4a35-9d5f-0f734e46ef03
# ╠═e33fe5dd-db88-4de9-8816-b4c36d6765d1
# ╟─b2e98c3f-5de5-4c23-bbb7-80d39cf219fd
# ╠═fb6f9b74-1687-4846-b070-58581a1b234b
# ╟─4bb856a1-5696-42c8-8e63-fd7abe0eb276
# ╠═9fa86072-9483-4b63-8dc9-e94469c1a8e9
# ╠═3367e5bd-3159-4327-90a4-01b7799e7064
# ╠═c30df4ad-de11-446c-84ed-33755a78e5aa
# ╠═e7284fe8-79b8-4684-922c-f004cb98333e
