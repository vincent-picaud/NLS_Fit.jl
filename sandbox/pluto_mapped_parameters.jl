### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 2a9bf136-9e12-11ec-2a2b-b7fdaa69775a
begin
	import Pkg
	
	using PyPlot
	using Random
	
	# Pkg.develop(path="/home/picaud/GitHub/NLS_Fit.jl")
	Pkg.develop(path="..")
	using NLS_Fit
end

# ╔═╡ 23675958-4f1b-4357-ae41-13613579414d
md"""
# Mapped parameters

This demo shows to define 3 peaks sharing an affine law for their shape factor σ
"""

# ╔═╡ 00481ef5-e39d-418e-8e0c-acf2b86691d6
md"""
## Synthetic data
"""

# ╔═╡ a7bfe845-5798-4286-a0ce-b2c6c0197809
md"""
The model is the sum of 3 Gaussain peaks. Parameters are as follows:
```math
\begin{array}{l ccc l}
\theta = [ & h_1 & \mu_1 & \sigma_1 &  h_2 & \mu_2 & \sigma_2 & h_3 & \mu_3 & \sigma_3 & ]\\
 & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9
\end{array}
```
"""

# ╔═╡ d97b029d-214b-440d-a952-cdf46076c98a
begin
	Random.seed!(1234)

	model = Gaussian_Peak() + Gaussian_Peak() + Gaussian_Peak()

	h₁, μ₁, σ₁ = 1, 5, 1
	h₂, μ₂, σ₂ = 1, 10, 1.5
	h₃, μ₃, σ₃ = 1, 20, 2
	
	θ = Float64[ h₁, μ₁, σ₁, h₂, μ₂, σ₂, h₃, μ₃, σ₃]
end

# ╔═╡ dbc94276-7474-4d36-a4c7-66415e61e9c5
begin
	X=Float64[1:0.25:30;]
	n=length(X)
	Y=eval_y(model,X,θ) + 0.1*(rand(n) .- 0.5)

	fig, ax = PyPlot.subplots()
	ax.plot(X,Y, "+", label="Raw")
	fig
end

# ╔═╡ 8064dd97-d036-4cc4-9bec-d21f497d487e
md"""
# Define the model to fit
"""

# ╔═╡ 1542f565-be72-48db-969f-bcbb5f8f4aa9
md"""
## Define σ law
"""

# ╔═╡ 4e7ac74a-b12b-4718-9c01-4597c8cc1836
md"""
The idea is to transform the initial θ parameter vector:

```math
\begin{array}{lcccccccccl}
\theta = [ & h_1 & \mu_1 & \sigma_1 &  h_2 & \mu_2 & \sigma_2 & h_3 & \mu_3 & \sigma_3 & ]\\
 & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9
\end{array}
```

into

```math
\hat{\theta} = 
\begin{array}{l ccccc l}
[ & h_1 & \mu_1 & h_2 & \mu_2 & h_3 & \mu_3 & ]
\end{array}
\textrm{ with }
\left\{
\begin{array}{lcr}
\sigma_1 & = & f(\mu_1) \\
\sigma_2 & = & f(\mu_2) \\
\sigma_3 & = & f(\mu_3) \\
\end{array}
\right.
```
Thus, we apply $f$ from components $2,5,8$ to components $3,6,9$. 

This operation is performed thanks to the `create_model_transform_src_insert_dest` function:
"""

# ╔═╡ 6cf3d077-a15c-43ca-b8fa-dce4ed5a6484
model_σ_μ, θ̂ = let
	f_σ_μ = Map_Affine(1.0  => 1.0, 30.0 => 5.0)
	src=[2,5,8]
	dest=[3,6,9]
	
	model_σ_μ = create_model_transform_src_insert_dest(model,f_σ_μ,src=>dest)
	θ̂ = Float64[h₁, μ₁, h₂, μ₂, h₃, μ₃, 1,1]

	model_σ_μ, θ̂
end

# ╔═╡ ea73741d-bdd7-474e-9c90-df609dede9b1


# ╔═╡ 03e7f62d-57ea-4c63-8fd8-ff51c6e8201d
begin
	Ŷ=eval_y(model_σ_μ, X, θ̂)
	ax.plot(X,Ŷ,label="mlkm")
	fig
end

# ╔═╡ c5ad04ca-cd6d-46d3-8d8f-9added482723
md"""
One can now modify $\hat{\theta}$ to impose constant position $\mu$.

From 
```math
\begin{array}{lccccccl}
\hat{\theta} = [ & h_1 & \mu_1 & h_2 & \mu_2 & h_3 & \mu_3 & ]\\
 & 1 & 2 & 3 & 4 & 5 & 6 &
\end{array}
```
we define 
```math
\begin{array}{lcccl}
\tilde{\theta} = [ & h_1 & h_2 & h_3 & ]
\end{array}
\textrm{ with }
\left\{
\begin{array}{lcr}
\mu_1 & = & 5 \\
\mu_2 & = & 10 \\
\mu_3 & = & 20 \\
\end{array}
\right.
```
"""

# ╔═╡ 23225877-308a-4cd5-9308-00d80ae43f94
model_h, θ̃ = let
	indices = [2, 4, 6]
	values  = [5, 10, 20]
	
	model_h = Model2Fit_Const_Parameters(model_σ_μ,indices,values)
	θ̃ = Float64[h₁, h₂, h₃, 1, 1]

	model_h, θ̃
end

# ╔═╡ 533e5948-65ea-4b25-9b6f-c8ce26213022
md"""
## Solve the problem
"""

# ╔═╡ 7c856e7c-cc43-4e9e-9261-49160e937946
result = let
	nls = NLS_ForwardDiff_From_Model2Fit(model_h,X,Y)
	conf = NLS_Solver.LevenbergMarquardt_Conf()
	result = NLS_Solver.solve(nls,θ̃,conf)
end

# ╔═╡ ff21fd91-0afc-4a66-abea-fdeaeab41a95
θ_solution = NLS_Solver.solution(result)

# ╔═╡ 927c5dd4-8ce0-47da-9ec9-c26c37331619
md"""
## Plot solution
"""

# ╔═╡ 34510940-986a-42e7-a119-7748493bab94
begin
	Y_solution=eval_y(model_h, X, θ_solution)
	ax.plot(X,Y_solution,label="mlkm")
	fig
end

# ╔═╡ 2ec887e4-42a7-4944-9906-5468a6077b57
md"""
TODO: retrieve initial σ
"""

# ╔═╡ Cell order:
# ╠═2a9bf136-9e12-11ec-2a2b-b7fdaa69775a
# ╠═23675958-4f1b-4357-ae41-13613579414d
# ╠═00481ef5-e39d-418e-8e0c-acf2b86691d6
# ╟─a7bfe845-5798-4286-a0ce-b2c6c0197809
# ╟─d97b029d-214b-440d-a952-cdf46076c98a
# ╟─dbc94276-7474-4d36-a4c7-66415e61e9c5
# ╟─8064dd97-d036-4cc4-9bec-d21f497d487e
# ╟─1542f565-be72-48db-969f-bcbb5f8f4aa9
# ╟─4e7ac74a-b12b-4718-9c01-4597c8cc1836
# ╠═6cf3d077-a15c-43ca-b8fa-dce4ed5a6484
# ╠═ea73741d-bdd7-474e-9c90-df609dede9b1
# ╠═03e7f62d-57ea-4c63-8fd8-ff51c6e8201d
# ╟─c5ad04ca-cd6d-46d3-8d8f-9added482723
# ╠═23225877-308a-4cd5-9308-00d80ae43f94
# ╠═533e5948-65ea-4b25-9b6f-c8ce26213022
# ╠═7c856e7c-cc43-4e9e-9261-49160e937946
# ╠═ff21fd91-0afc-4a66-abea-fdeaeab41a95
# ╟─927c5dd4-8ce0-47da-9ec9-c26c37331619
# ╠═34510940-986a-42e7-a119-7748493bab94
# ╠═2ec887e4-42a7-4944-9906-5468a6077b57
