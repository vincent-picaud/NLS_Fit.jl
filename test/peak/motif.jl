@testset "motif.jl" begin
    n = 5

    profile=rand(n,2)
    
    model = Peak_Motif(Gaussian_Peak(),profile)
    
    @test parameter_size(model) == 2

    θ = rand(parameter_size(model))

    @test (@benchmark eval_y($model,0.0,$θ)).allocs == 0
end
