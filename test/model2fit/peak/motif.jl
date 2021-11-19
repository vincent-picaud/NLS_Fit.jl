@testset "motif.jl" begin
    n = 5

    profile=rand(n,2)
    
    model = Peak_Motif(Gaussian_Peak(),profile)
    
    @test parameter_size(model) == 2

    θ = rand(parameter_size(model))

    X = rand(10)
    Y = zeros(10)
    
    @test (@benchmark eval_y!($model,$Y,$X,$θ)).allocs == 0
end
