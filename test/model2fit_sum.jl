@testset "model2fit_sum.jl" begin

    @testset "Empty left" begin
        
        model = Model2Fit_Empty()
        model = model + Gaussian_Peak()

        @test parameter_size(model) == 3

        θ=Float64[1, 0, 1]

        @test eval_y(model,[0.0],θ) ≈ [1.0]

        X = rand(10)
        Y = alloc_y(model,X,θ)
        
        @test (@benchmark eval_y!($model, $Y, $X,$θ)).allocs == 0
    end


    @testset "Empty right" begin
        
        model = Model2Fit_Empty()
        model = Gaussian_Peak() + model

        @test parameter_size(model) == 3

        θ=Float64[1, 0, 1]

        @test eval_y(model,[0.0],θ) ≈ [1.0]

        X = rand(10)
        Y = alloc_y(model,X,θ)
        
        @test (@benchmark eval_y!($model, $Y, $X,$θ)).allocs == 0
        
    end

    @testset "2 terms sum" begin
        
        model = Gaussian_Peak() + Gaussian_Peak()

        @test parameter_size(model) == 6

        θ=Float64[1, 0, 1, 1, 0, 4]

        @test eval_y(model,[0.0],θ) ≈ [2.0]

        X = rand(10)
        Y = alloc_y(model,X,θ)
        
        @test (@benchmark eval_y!($model, $Y, $X,$θ)).allocs == 0
        
    end

    @testset "3 terms sum" begin
        
        model = Gaussian_Peak() + Gaussian_Peak() + Gaussian_Peak()

        @test parameter_size(model) == 9

        θ=Float64[1, 0, 1, 1, 0, 4, 1, 0, 4]

        @test eval_y(model,[0.0],θ) ≈ [3.0]

        X = rand(10)
        Y = alloc_y(model,X,θ)
        
        @test (@benchmark eval_y!($model, $Y, $X,$θ)).allocs == 0
        
    end

end
