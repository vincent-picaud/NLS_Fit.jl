@testset "model2fit_sum.jl" begin

    @testset "Sum of 2 Gaussian with const centers" begin

        model = Gaussian_Peak() + Gaussian_Peak()
        model_θ = Float64[1,5,2,  # h1=1,μ1=5,σ1=2
                          2,10,3] # h2=2,μ2=10,σ2=3

        X = collect(1:1.0:10)
        Y = eval_y(model,X,model_θ)

        indices_μ = [2,5]
        values_μ = [5.0,10.0]
        model_const_μ = Model2Fit_Const_Parameters(model,[2,5],[5.0,10.0])
        model_const_μ_θ = deleteat!(copy(model_θ),indices_μ)
        
        @test parameter_size(model_const_μ) == length(model_const_μ_θ)

        @test get_model(model_const_μ) == model
        @test get_model_θ(model_const_μ, model_const_μ_θ) == model_θ
        
        Y_const = eval_y(model_const_μ,X,model_const_μ_θ)
        @test Y == Y_const

    end

end
