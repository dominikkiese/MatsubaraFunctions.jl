module MatsubaraFunctions

    using PrecompileTools
    @recompile_invalidations begin
        using LinearAlgebra
        using StaticArrays
        using MPI
        using Polyester
        using HDF5
        using Aqua
        using Documenter
        using Coverage
    end

    # macro to enable debug mode
    DEBUG() = false

    macro DEBUG(expr, msgs)
        esc(:(if $(@__MODULE__).DEBUG() @assert($expr, $msgs...) end))
    end

    # dummy functions, use if not overloaded
    index(p)             = p
    mesh_index(p, m)     = index(p)
    mesh_index_bc(p, m)  = mesh_index(p, m)
    is_inbounds_bc(p, m) = true

    include("mesh/mesh.jl")
    include("func/func.jl")
    include("misc/mpi_helpers.jl")
    include("misc/pade.jl")
    include("misc/triqs_interface.jl")
end