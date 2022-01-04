# TODO: add tests


# An helper private function
#
# Role: modify src index taken into account that dest component are removed
#
# Algo: for each dest index inferior to a given src index, src index
# must be reduced by one unit.
#
function src_indices_after_dest_removed(src::AbstractVector{Int},
                                        dest::AbstractVector{Int})
    @assert issorted(dest) # only dest has to be sorted (reason: the
                           # insertion algo requires an ordered dest)
    @assert length(src)==length(dest)
    
    n = length(src)
    
    count =  zeros(Int,n)
    i_count = 1

    # note this works even if src, dest are not sorted but this has
    # n^2 complexity
    for src_i in src
        for dest_i in dest
            if dest_i < src_i
                count[i_count] += 1
            else
                @assert dest_i != src_i "this makes no sense"
            end
        end 
        i_count += 1
    end
    src - count
end

@doc raw"""

Define a map ``g`` of the form
```math
\theta = g_{\hat{\theta}_f}(\hat{\theta}_m)
```

Where the resulting ``\theta`` vector is first initialized to
``\hat{\theta}_m``. This vector is then completed by insertions: for
each ``i\in\text{src}`` one computes

```math
f_{\hat{\theta}_f}(`\theta[i])
```
and inserts the result at position `j=dest[i]`.

# Constructor

```julia
Transformed_Parameter_Src_Dest_Map(f_map, src=>dest)
```
- `f_map` is the stored [`Abstract_Map`](@ref) map,

- `src=>dest` is a pair of two `Int` vectors used to define
  indices. These indices are assumed to be sorted. Also, by
  construction, we must have `src ∩ dest = ∅`. (the reason is that as
  dest indices are removed, there would no src component to apply
  ``f``)

# Example

TODO: but before add separable map 
"""
struct Transformed_Parameter_Src_Dest_Map <: Abstract_Map
    _f_map::Abstract_Map
    _src_dest::Matrix{Int} # src=col(1), dest=col(2)
end

function Transformed_Parameter_Src_Dest_Map(f_map::Abstract_Map, src_dest::Pair{SOURCE,DEST}) where {SOURCE<:AbstractVector{Int},DEST<:AbstractVector{Int}}

    (src, dest) = src_dest

    @assert length(src) == length(dest)

    # Form [src,dest] matrix & sort rows wrt dest
    src_dest = hcat(src, dest)
    src_dest = sortslices(src_dest, dims = 1, by = x -> last(x))

    # Update src indices
    src_dest[:, 1] .= src_indices_after_dest_removed(src_dest[:, 1], src_dest[:, 2])

    Transformed_Parameter_Src_Dest_Map(f_map, src_dest)
end

parameter_size(tm::Transformed_Parameter_Src_Dest_Map) = parameter_size(tm._f_map)

function eval_map(tm::Transformed_Parameter_Src_Dest_Map,hat_model_θ::AbstractVector,hat_map_θ::AbstractVector)
    @assert parameter_size(tm) == length(hat_map_θ)

    src  = @view tm._src_dest[:,1]
    dest = @view tm._src_dest[:,2]
    
    θf = eval_map(tm._f_map,(@view hat_model_θ[src]),hat_map_θ)
   
    θ = insert_some_elements(hat_model_θ, dest, θf)

    θ
end
