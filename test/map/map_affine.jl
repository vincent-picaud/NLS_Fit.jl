@testset "affine.jl" begin

    map = NLS_Fit.Map_Affine(1=>2,3=>5)

    @test parameter_size(map) == 2

    θ = Float32[1, 1]
    X_hat = Int[1,2,3]
    X = eval_map(map,X_hat,θ)

    @test eltype(X) == Float32
    @test X ≈ Float32[2,3.5,5]
    
end

