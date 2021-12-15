# insert some elements: useful when one wants to modify θ by inserting
# new parameters
#

export insert_some_elements

# ****************************************************************

#
# How to do that with issorted is not clear:
#
# https://discourse.julialang.org/t/understanding-issorteds-lt-keyword/60918/42?page=3
#
# ... to avoid any missuses of issorted we define our own test
#
# This function returns true in case of _strict_ ordering, false otherwise
#
function is_strictly_ordered(indices::AbstractVector{Int})

    iter = iterate(indices)
    if iter == nothing return true end # empty sequence case

    iter_next = iterate(indices,last(iter))

    while iter_next != nothing
        first(iter) >= first(iter_next) && return false
        iter=iter_next
        iter_next=iterate(indices,last(iter_next))
    end

    true 
end


# A priori, can easily be converted for arbitrary axes:
# - remplace i_dest by iterate(X_dest)
# - remplace i_src by iterate(X)
# then
# - remplace i_dest += 1 by iterate(X_dest,last(i_dest))
# - same idea for i_src += 1
#
@doc raw"""
```julia
insert_some_elements(X::AbstractVector{T},
                     indices::AbstractVector{Int},
                     elements::AbstractVector) -> AbstractVector{T}
```

Insert `elements` in `X` at positions `indices`.

Note: `indices` must be strictly ordered.

Example: this functions does a work inverse to the `deleteat!`

```jldoctest
julia> X = Any[1:10;]
10-element Vector{Any}:
  1
  2
  3
  4
  5
  6
  7
  8
  9
 10

julia> indices = Int[2, 3, 5, 6, 9, 10];

julia> elements = Float64[2, 3, 5, 6, 9, 10];

julia> deleteat!(X, indices)
4-element Vector{Any}:
 1
 4
 7
 8

julia> insert_some_elements(X,indices,elements)
10-element Vector{Any}:
  1
  2.0
  3.0
  4
  5.0
  6.0
  7
  8
  9.0
 10.0

```
"""
function insert_some_elements end

# We first define a helper where elements type is free
# This allows to use the methods with *iterable* elements
#
function _insert_some_elements(X::AbstractVector{T},
                              indices::AbstractVector{Int},
                              elements) where {T} # elements must be iterable

    @assert is_strictly_ordered(indices)
    @assert length(indices)==0 || 1 ≤ minimum(indices) ≤ length(X)+length(indices)

    X_dest = Vector{T}(undef,length(X)+length(indices))

    i_dest = 1
    i_src = 1
    for (i,e) in zip(indices,elements)
        while i_dest<i
            X_dest[i_dest]=X[i_src]
            i_dest += 1
            i_src += 1
        end
        X_dest[i_dest]=e
        i_dest += 1
    end
    
    n = length(X)
    while i_src<=n
        X_dest[i_dest]=X[i_src]
        i_dest += 1
        i_src += 1
    end
  
    X_dest
end

# Specialize when elements is an AbstractVector
#
function insert_some_elements(X::AbstractVector{T},
                              indices::AbstractVector{Int},
                              elements::AbstractVector) where {T}

    @assert all(V->axes(V,1) isa Base.OneTo,(X,indices,elements)) # V[i], with i=1..n
    @assert length(indices)==length(elements)

    _insert_some_elements(X,
                          indices,
                          elements)
end

# Specialize when elements is a constant
#
function insert_some_elements(X::AbstractVector{T},
                              indices::AbstractVector{Int},
                              elements::T) where {T}

    _insert_some_elements(X,
                          indices,
                          Base.Iterators.repeated(elements))
end

