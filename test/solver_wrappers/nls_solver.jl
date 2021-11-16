@testset "nls_solver.jl" begin

    using NLS_Solver
    
    @testset "memory allocs" begin
 
        n = 10
        X = Float64[1:n;]
        Y = rand(n)
        
        model = Gaussian_Peak()
        
        θ = Float64[1,0,1]
        
        nls = NLS_ForwardDiff_From_Model2Fit(model,X,Y)

        @test NLS_Solver.parameter_size(nls) == 3
        
    end
    
end
