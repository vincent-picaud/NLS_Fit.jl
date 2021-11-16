@testset "model2fit_sum.jl - empty left" begin
    
    model = Model2Fit_Empty()
    model = model + Gaussian_Peak()

    @test parameter_size(model) == 3

    θ=Float64[1, 0, 1]

    @test eval_y(model,0.0,θ) ≈ 1

    @test (@allocated eval_y(model,0.0,θ)) == 0
    
end

@testset "model2fit_sum.jl - empty right" begin
    
    model = Model2Fit_Empty()
    model = Gaussian_Peak() + model

    @test parameter_size(model) == 3

    θ=Float64[1, 0, 1]

    @test eval_y(model,0.0,θ) ≈ 1

    @test (@allocated eval_y(model,0.0,θ)) == 0
    
end

@testset "model2fit_sum.jl - 2" begin
    
    model = Gaussian_Peak() + Gaussian_Peak()

    @test parameter_size(model) == 6

    θ=Float64[1, 0, 1, 1, 0, 4]

    @test eval_y(model,0.0,θ) ≈ 2

    @test (@allocated eval_y(model,0.0,θ)) == 0
    
end

@testset "model2fit_sum.jl - 3" begin
    
    model = Gaussian_Peak() + Gaussian_Peak() + Gaussian_Peak()

    @test parameter_size(model) == 9

    θ=Float64[1, 0, 1, 1, 0, 4, 1, 0, 4]

    @test eval_y(model,0.0,θ) ≈ 3

    @test (@allocated eval_y(model,0.0,θ)) == 0
    
end
