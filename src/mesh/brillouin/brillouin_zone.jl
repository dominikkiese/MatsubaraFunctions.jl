# type def and accessors
#-------------------------------------------------------------------------------#

"""
    struct BrillouinZone{N, P}

BrillouinZone struct with fields:
* `L         :: Int`                       : linear system size
* `basis     :: SMatrix{N, N, Float64, P}` : matrix with reciprocal lattice vectors as columns
* `inv_basis :: SMatrix{N, N, Float64, P}` : matrix inverse of `basis`
"""
struct BrillouinZone{N, P}
    L         :: Int
    basis     :: SMatrix{N, N, Float64, P}
    inv_basis :: SMatrix{N, N, Float64, P}

    function BrillouinZone(L :: Int, basis :: SMatrix{N, N, Float64, P}) where {N, P}
        return new{N, P}(L, basis, inv(basis))
    end 

    function BrillouinZone(L :: Int, vecs :: Vararg{SVector{N, Float64}, N}) where {N}
        return BrillouinZone(L, hcat(vecs...))
    end 
end

"""
    function basis(bz :: BrillouinZone{N, P}) :: SMatrix{N, N, Float64, P} where {N, P} 

Returns `bz.basis`
"""
function basis(bz :: BrillouinZone{N, P}) :: SMatrix{N, N, Float64, P} where {N, P} 
    return bz.basis 
end

"""
    function basis(bz  :: BrillouinZone{N, P}, idx :: Int) :: SVector{N, Float64} where {N, P}

Returns basis vector with index `idx`
"""
function basis(bz  :: BrillouinZone{N, P}, idx :: Int) :: SVector{N, Float64} where {N, P}
    return bz.basis[:, idx]
end

"""
    function inv_basis(bz :: BrillouinZone{N, P}) :: SMatrix{N, N, Float64, P} where {N, P} 

Returns `bz.inv_basis`
"""
function inv_basis(bz :: BrillouinZone{N, P}) :: SMatrix{N, N, Float64, P} where {N, P} 
    return bz.inv_basis 
end

"""
    function inv_basis(bz  :: BrillouinZone{N, P}, idx :: Int) :: SVector{N, Float64} where {N, P}

Returns inverse basis vector with index `idx`
"""
function inv_basis(bz  :: BrillouinZone{N, P}, idx :: Int) :: SVector{N, Float64} where {N, P}
    return bz.inv_basis[:, idx]
end

# conversion from reciprocal to euclidean coordinates
#-------------------------------------------------------------------------------#

"""
    function euclidean(k :: BrillouinPoint{N}, bz :: BrillouinZone{N, P}) :: SVector{N, Float64} where {N, P}

Convert reciprocal to euclidean coordinates
"""
function euclidean(k :: BrillouinPoint{N}, bz :: BrillouinZone{N, P}) :: SVector{N, Float64} where {N, P}
    return basis(bz) * (value(k) ./ bz.L)
end

# conversion from euclidean to reciprocal coordinates
#-------------------------------------------------------------------------------#

"""
    function reciprocal(k :: T, bz :: BrillouinZone{N, P}) :: SVector{N, Float64} where {N, P, T <: AbstractVector{Float64}}

Convert euclidean to reciprocal coordinates
"""
function reciprocal(k :: T, bz :: BrillouinZone{N, P}) :: SVector{N, Float64} where {N, P, T <: AbstractVector{Float64}}
    @DEBUG length(k) == N "Length mismatch for input vector"
    return inv_basis(bz) * (bz.L .* k)
end

# bounds checking
#-------------------------------------------------------------------------------#

"""
    function is_inbounds(k :: BrillouinPoint{N}, bz :: BrillouinZone{N, P}) :: Bool where {N, P}

Checks if reciprocal coordinates in bounds
"""
function is_inbounds(k :: BrillouinPoint{N}, bz :: BrillouinZone{N, P}) :: Bool where {N, P}
    return all(kn -> 0 <= kn < bz.L, value(k))
end

"""
    function is_inbounds(k :: T, bz :: BrillouinZone{N, P}) :: Bool where {N, P, T <: AbstractVector{Float64}}

Checks if euclidean coordinates in bounds
"""
function is_inbounds(k :: T, bz :: BrillouinZone{N, P}) :: Bool where {N, P, T <: AbstractVector{Float64}}
    @DEBUG length(k) == N "Length mismatch for input vector"
    return all(kn -> 0 <= kn < bz.L, reciprocal(k, bz))
end

# periodic boundary conditions
#-------------------------------------------------------------------------------#

function positive_modulo(idx :: T, L; thresh = T(0)) :: T where {T}
    idxp = T(idx % L) 
    return idxp >= thresh ? idxp : idxp + L 
end

"""
    function fold_back(k :: BrillouinPoint{N}, bz :: BrillouinZone{N, P}) :: BrillouinPoint{N} where {N, P}

Use periodic boundary conditions to fold `k` back into mesh
"""
function fold_back(k :: BrillouinPoint{N}, bz :: BrillouinZone{N, P}) :: BrillouinPoint{N} where {N, P}
    return BrillouinPoint(ntuple(n -> positive_modulo(k[n], bz.L), N)...)
end

"""
    function fold_back(k :: T, bz :: BrillouinZone{N, P}) :: SVector{N, Float64} where {N, P, T <: AbstractVector{Float64}}

Use periodic boundary conditions to fold `k` back into mesh
"""
function fold_back(k :: T, bz :: BrillouinZone{N, P}) :: SVector{N, Float64} where {N, P, T <: AbstractVector{Float64}}
    # slower than back folding in reciprocal space, should be avoided if possible
    @DEBUG length(k) == N "Length mismatch for input vector"

    x = reciprocal(k, bz)
    return basis(bz) * (SVector{N, Float64}(ntuple(n -> positive_modulo(x[n], bz.L; thresh = -1e-14), N)...) ./ bz.L)
end

# mapping to Wigner-Seitz cell
#-------------------------------------------------------------------------------#

"""
    function get_shifts(bz :: BrillouinZone{N, P}) :: Matrix{SVector{N, Float64}} where {N, P}

Generate all primitive shifts in euclidean space from Brillouin zone
"""
function get_shifts(bz :: BrillouinZone{N, P}) :: Matrix{SVector{N, Float64}} where {N, P}
    iters = Iterators.product(ntuple(n -> -1 : 1, N)...)
    return SVector{N, Float64}[basis(bz) * SVector{N, Float64}(iter...) for iter in iters]
end

"""
    function to_Wigner_Seitz(k :: T, bz :: BrillouinZone{N, P}) :: SVector{N, Float64} where {N, P, T <: AbstractVector{Float64}}

Map k to Wigner Seitz cell at Γ = 0
"""
function to_Wigner_Seitz(k :: T, bz :: BrillouinZone{N, P}) :: SVector{N, Float64} where {N, P, T <: AbstractVector{Float64}}
    @DEBUG length(k) == N "Length mismatch for input vector"

    reducible = true 
    shifts    = get_shifts(bz)
    kp        = fold_back(k, bz) # start at k and fold back into primitive zone

    while reducible 
        reducible = false
        
        for shift in shifts
            kpp = kp - shift

            if norm(kp) - norm(kpp) > 1e-14
                kp        = kpp 
                reducible = true
            end 
        end 
    end 

    return kp 
end 

"""
    function to_Wigner_Seitz(k :: T, shifts :: Matrix{SVector{N, Float64}}) :: SVector{N, Float64} where {N, T <: AbstractVector{Float64}}

Map k to Wigner Seitz cell at Γ = 0 using precomputed set of translations. Useful when mapping many points from the primitive zone.
"""
function to_Wigner_Seitz(k :: T, shifts :: Matrix{SVector{N, Float64}}) :: SVector{N, Float64} where {N, T <: AbstractVector{Float64}}
    @DEBUG length(k) == N "Length mismatch for input vector"

    reducible = true 
    kp        = k 

    while reducible 
        reducible = false
        
        for shift in shifts
            kpp = kp - shift

            if norm(kp) - norm(kpp) > 1e-14
                kp        = kpp 
                reducible = true
            end 
        end 
    end 

    return kp 
end

# info
#-------------------------------------------------------------------------------#

function info(bz :: BrillouinZone{N, P}) where {N, P}
    println(CYAN, BOLD, "BrillouinZone ", RESET, "of dimension ", CYAN, BOLD, "$(N)", RESET)
    println("=> L     : $(bz.L)")
    println("=> basis : $(basis(bz))")
    return nothing 
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
    get_shifts,
    to_Wigner_Seitz,
    info