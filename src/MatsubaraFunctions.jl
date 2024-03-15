module MatsubaraFunctions

    using PrecompileTools

    @recompile_invalidations begin
        using Printf
        using LinearAlgebra
        using StaticArrays
        using MPI
        using HDF5
        using Aqua
        using Documenter
    end

    # macro to unlock debug mode
    DEBUG() = false

    macro DEBUG(expr, msgs)
        esc(:(if $(@__MODULE__).DEBUG() @assert($expr, $msgs...) end))
    end

    # dummy functions, use if not overloaded
    mesh_index(p, m)     = p
    mesh_index_bc(p, m)  = mesh_index(p, m)
    is_inbounds_bc(p, m) = true

    include("mesh/mesh.jl")
    include("func/func.jl")

    include("misc/mpi_helpers.jl")
    include("misc/pade.jl")
    include("misc/pulay.jl")
end