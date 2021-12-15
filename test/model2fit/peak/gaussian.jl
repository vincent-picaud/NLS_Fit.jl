@testset "gaussian.jl" begin
    model = Gaussian_Peak()

    @test parameter_size(model) == 3

    θ=Float64[1, 0, 1]

    @test eval_y(model,[0.0],θ) ≈ [1.0]

    X=rand(10)
    Y=alloc_y(model,X,θ)
    
    @test (@benchmark accumulate_y!($model,$Y,$X,$θ)).allocs == 0
    
end
