@testset "insert_some_elements.jl" begin
    @testset "elements = vector" begin
        using Random: seed!
        using StatsBase: sample
        
        seed!(1234)

        for i in 0:10
            X=[1:10;]
            indices = sort(sample(1:length(X),i,replace=false))
            elements=X[indices]
            
            X_deflated = deleteat!(copy(X),indices)
            X_recons = insert_some_elements(X_deflated,indices,elements)

            @test X_recons == X
        end
    end

    @testset "elements = constant" begin
        using Random: seed!
        using StatsBase: sample
        
        seed!(1234)

        for i in 0:10
            X=[1:10;]
            indices = sort(sample(1:length(X),i,replace=false))
            element = 2021
            X[indices] .= element
            
            X_deflated = deleteat!(copy(X),indices)
            X_recons = insert_some_elements(X_deflated,indices,element)

            @test X_recons == X
        end
    end
end
