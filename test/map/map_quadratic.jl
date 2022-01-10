@testset "quadratic.jl" begin

    @testset "scaled" begin

        p(x)=1+2*x+3*x*x
        
        map = NLS_Fit.Map_Quadratic(1=>p(1),3=>p(3),10=>p(10))

        @test parameter_size(map) == 3

        θ = Float32[1, 1, 1]
        X_hat = Int[1,2,3,4,10]
        X = eval_map(map,X_hat,θ)

        @test eltype(X) == Float32
        @test X ≈ p.(X_hat)
        
    end

    @testset "unscaled" begin
        
        p(x) = 3+2*x+1*x*x
        
        map = NLS_Fit.Map_Quadratic(1,3,10)
        
        @test parameter_size(map) == 3
        
        θ = Float32[p(1), p(3), p(10)]
        X_hat = Int[1,2,3,4,10]
        X = eval_map(map,X_hat,θ)

        @test eltype(X) == Float32
        @test X ≈ p.(X_hat)
        
    end
end

