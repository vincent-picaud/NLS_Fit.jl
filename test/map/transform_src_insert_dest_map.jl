@testset "transform_src_insert_dest_map.jl" begin

    @testset "test 1" begin

        f(X,θ̂) = θ̂[1] .+ θ̂[2]*X

        f_map = Map_From_VectFunc(2,f)
        hat_θ_map = Float64[0,-1] # change sign

        # insert -θ4, at position 2
        # insert -θ1, at position 3
        src  = [4,1]
        dest = [2,3]

        # delete insert dest [1,2,3,4,5] -> [1,4,5]
        hat_θ_model = Float64[1:5;]
        deleteat!(hat_θ_model,dest)

        map_s_d = NLS_Fit.Transform_Src_Insert_Dest_Map(f_map,src=>dest)

        θ = eval_map(map_s_d,hat_θ_model,hat_θ_map)

        @test parameter_size( map_s_d) == 2
        @test θ ≈ Float64[1,-4,-1,4,5]
    end
    
    @testset "test 2" begin
        # same as test 1 but with permuted dest
        
        f(X,θ̂) = θ̂[1] .+ θ̂[2]*X
        
        f_map = Map_From_VectFunc(2,f)
        hat_θ_map = Float64[0,-1] # change sign
        
        src  = [4,1]
        dest = [3,2]
        
        hat_θ_model = Float64[1:5;]
        deleteat!(hat_θ_model,sort(dest))
        
        map_s_d = NLS_Fit.Transform_Src_Insert_Dest_Map(f_map,src=>dest)
        
        θ = eval_map(map_s_d,hat_θ_model,hat_θ_map)
        
        @test parameter_size( map_s_d) == 2
        @test θ ≈ Float64[1,-1,-4,4,5]
    end 
end
