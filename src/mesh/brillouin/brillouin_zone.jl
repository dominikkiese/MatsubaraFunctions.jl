# type def and accessors
#-------------------------------------------------------------------------------#

"""
    struct BrillouinZone{N}

BrillouinZone struct with fields:
* `L         :: Int64`                  : linear system size
* `basis     :: SMatrix{N, N, Float64}` : matrix with reciprocal lattice vectors as columns
* `inv_basis :: SMatrix{N, N, Float64}` : inverse of `basis`
"""
struct BrillouinZone{N}
    L         :: Int64
    basis     :: SMatrix{N, N, Float64}
    inv_basis :: SMatrix{N, N, Float64}

    function BrillouinZone(
        L     :: Int64,
        basis :: SMatrix{N, N, Float64}
        )     :: BrillouinZone{N} where {N}

        return new{N}(L, basis, inv(basis))
    end 

    function BrillouinZone(
        L    :: Int64,
        vecs :: Vararg{SVector{N, Float64}, N}
        )    :: BrillouinZone{N} where {N}

        return BrillouinZone(L, hcat(vecs...))
    end 
end

"""
    function basis(
        bz :: BrillouinZone{N}
        )  :: SMatrix{N, N, Float64} where {N} 

Returns `bz.basis`
"""
function basis(
    bz :: BrillouinZone{N}
    )  :: SMatrix{N, N, Float64} where {N}

    return bz.basis 
end

"""
    function basis(
        bz  :: BrillouinZone{N},
        idx :: Int64
        )   :: SVector{N, Float64} where {N}

Returns basis vector with index `idx`
"""
function basis(
    bz  :: BrillouinZone{N},
    idx :: Int64
    )   :: SVector{N, Float64} where {N}

    return bz.basis[:, idx]
end

"""
    function inv_basis(
        bz :: BrillouinZone{N}
        )  :: SMatrix{N, N, Float64} where {N} 

Returns `bz.inv_basis`
"""
function inv_basis(
    bz :: BrillouinZone{N}
    )  :: SMatrix{N, N, Float64} where {N} 

    return bz.inv_basis 
end

"""
    function inv_basis(
        bz  :: BrillouinZone{N},
        idx :: Int64
        )   :: SVector{N, Float64} where {N}

Returns inverse basis vector with index `idx`
"""
function inv_basis(
    bz  :: BrillouinZone{N},
    idx :: Int64
    )   :: SVector{N, Float64} where {N}

    return bz.inv_basis[:, idx]
end

# conversion from reciprocal to euclidean coordinates
#-------------------------------------------------------------------------------#

"""
    function euclidean(
        k  :: BrillouinPoint{N},
        bz :: BrillouinZone{N}
        )  :: SVector{N, Float64} where {N}

Convert reciprocal to euclidean coordinates
"""
function euclidean(
    k  :: BrillouinPoint{N},
    bz :: BrillouinZone{N}
    )  :: SVector{N, Float64} where {N}

    return basis(bz) * (index(k) ./ bz.L)
end

# conversion from euclidean to reciprocal coordinates
#-------------------------------------------------------------------------------#

"""
    function reciprocal(
        k  :: SVector{N, Float64},
        bz :: BrillouinZone{N}
        )  :: SVector{N, Float64} where {N}

Convert euclidean to reciprocal coordinates
"""
function reciprocal(
    k  :: SVector{N, Float64},
    bz :: BrillouinZone{N}
    )  :: SVector{N, Float64} where {N}

    return inv_basis(bz) * (bz.L .* k)
end

# bounds checking
#-------------------------------------------------------------------------------#

"""
    function is_inbounds(
        k  :: BrillouinPoint{N},
        bz :: BrillouinZone{N}
        )  :: Bool where {N}

Checks if reciprocal coordinates in bounds
"""
function is_inbounds(
    k  :: BrillouinPoint{N},
    bz :: BrillouinZone{N}
    )  :: Bool where {N}

    return all(kn -> 0 <= kn < bz.L, index(k))
end

"""
    function is_inbounds(
        k  :: SVector{N, Float64},
        bz :: BrillouinZone{N}
        )  :: Bool where {N}

Checks if euclidean coordinates in bounds
"""
function is_inbounds(
    k  :: SVector{N, Float64},
    bz :: BrillouinZone{N}
    )  :: Bool where {N}

    return all(kn -> 0 <= kn < bz.L, reciprocal(k, bz))
end

# periodic boundary conditions
#-------------------------------------------------------------------------------#

function positive_modulo(idx, L)
    idxp = idx % L 
    return idxp >= 0 ? idxp : idxp + L 
end

"""
    function fold_back(
        k  :: BrillouinPoint{N},
        bz :: BrillouinZone{N}
        )  :: BrillouinPoint{N} where {N}

Use periodic boundary conditions to fold `k` back into mesh
"""
function fold_back(
    k  :: BrillouinPoint{N},
    bz :: BrillouinZone{N}
    )  :: BrillouinPoint{N} where {N}

    return BrillouinPoint(map(kn -> positive_modulo(kn, bz.L), index(k)))
end

"""
    function fold_back(
        k  :: SVector{N, Float64},
        bz :: BrillouinZone{N}
        )  :: SVector{N, Float64} where {N}

Use periodic boundary conditions to fold `k` back into mesh
"""
function fold_back( # slower than back folding in reciprocal space, should be avoided if possible
    k  :: SVector{N, Float64},
    bz :: BrillouinZone{N}
    )  :: SVector{N, Float64} where {N}

    x = reciprocal(k, bz)
    return basis(bz) * (SVector{N, Float64}(ntuple(n -> positive_modulo(x[n], bz.L), N)...) ./ bz.L)
end

# mapping to Wigner-Seitz cell
#-------------------------------------------------------------------------------#

"""
    function to_Wigner_Seitz(
        k  :: SVector{N, Float64},
        bz :: BrillouinZone{N}
        )  :: SVector{N, Float64} where {N}

Map k to Wigner Seitz cell at Γ = 0
"""
function to_Wigner_Seitz(
    k  :: SVector{N, Float64},
    bz :: BrillouinZone{N}
    )  :: SVector{N, Float64} where {N}

    reducible = true 
    iters     = Iterators.product(ntuple(n -> -1 : 1, N)...)
    kp        = k 

    while reducible 
        reducible = false
        
        for iter in iters 
            kpp = kp - basis(bz) * SVector{N, Float64}(iter...)

            if norm(kp) - norm(kpp) > 1e-10
                kp        = kpp 
                reducible = true
            end 
        end 
    end 

    return kp 
end 

# export
#-------------------------------------------------------------------------------#

export 
    BrillouinZone,
    basis,
    inv_basis,
    euclidean,
    reciprocal,
    is_inbounds,
    fold_back,
    to_Wigner_Seitz