@testset "gaussian.jl" begin
    model = Gaussian_Peak()

    @test parameter_size(model) == 3

    θ=Float64[1, 0, 1]

    @test eval_y(model,0.0,θ) ≈ 1

    @test (@ballocated eval_y($model,0.0,$θ)) == 0
    
end
