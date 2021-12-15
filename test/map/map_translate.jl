@testset "translate.jl" begin

    @testset "scale=1" begin 

        map = Map_Translate()

        @test parameter_size(map) == 1

        θ = Float32[3]
        X_hat = Int[1,2,3]
        X = eval_map(map,X_hat,θ)

        @test eltype(X) == Float32
        @test X ≈ Float32[1,2,3] .+ (first(θ)-1)

    end

    @testset "scale=10" begin 

        map = Map_Translate(10.0)

        @test parameter_size(map) == 1

        θ = Float32[3]
        X_hat = Int[1,2,3]
        X = eval_map(map,X_hat,θ)

        @test eltype(X) == Float64
        @test X ≈ Float64[1,2,3] .+ 10 * (first(θ)-1)

    end
    
end
