@testset "transformed_parameters.jl" begin 

    @testset "Two Gaussians" begin

        # Two Gaussian peaks with σ(μ) an affine function
        #
        model = Gaussian_Peak() + Gaussian_Peak() 

        # model θ is:
        # idx:  1   2  3   4   5   6
        #  θ : h1, μ1, σ1, h2, μ2, σ2
        #
        # hence:
        #
        src  = [2,5]
        dest = [3,6]

        # The σ(μ) relation: σ(μ=1)=1 and σ(μ=100)=10 for [θ̂A, θ̂B]=[1, 1]
        #
        f_σ_μ = Map_Affine(1.0=>1.0,100.0=>10.0)

        # 
        g_map = NLS_Fit.Transformed_Parameter_Src_Dest_Map(f_σ_μ, src=>dest)

        model_σ_μ = Model2Fit_Transformed_Parameters(model,g_map)

        # Now the model_σ_μ parameters are:
        #
        # [ h1, μ1, h2, μ2, θ̂A, θ̂B ]
        #
        # the initial model is called with:
        #
        # [ h1, μ1, σ1 = affine(μ1, θ̂A, θ̂B), h2, μ2, , σ2 = affine(μ2, θ̂A, θ̂B) ]
        #
        model_σ_μ_θ = Float64[1,5, 2, 95, 1, 1]
        model_θ = get_model_θ(model_σ_μ, model_σ_μ_θ)

        @test model_θ ≈ [1.0, 5.0, 1.3636363636363638, 2.0, 95.0, 9.545454545454547]
        
        X = [1:10.0:101;]
        Y = eval_y(model_σ_μ, X, model_σ_μ_θ)

        @test Y ≈ [0.013538440136652369, 6.252150377485093e-5, 1.7807741296584817e-13, 3.4628155839082785e-10, 2.2470318744914768e-7, 4.865734936177968e-5, 0.003515984927922131, 0.08478224390932892, 0.6822163961931842, 1.8318872636421548, 1.6414765303163972]
        
    end

end 
