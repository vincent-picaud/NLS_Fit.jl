@testset "nls_solver.jl" begin

    @testset "memory allocs" begin

        n = 10
        X = Float64[1:n;]
        Y = rand(n)
        
        model = Gaussian_Peak()
        
        θ = Float64[1,0,1]
        
        nls = NLS_ForwardDiff_From_Model2Fit(model,X,Y)

        @test NLS_Solver.parameter_size(nls) == 3
        @test (@benchmark(NLS_Solver.eval_r($nls,$θ))).allocs == 1
        @test (@benchmark(NLS_Solver.eval_r_J($nls,$θ))).allocs == 9
        
    end
    
end
