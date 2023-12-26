module MatsubaraFunctions

    using PrecompileTools

    @recompile_invalidations begin
        using Printf
        using LinearAlgebra
        using StaticArrays
        using OffsetArrays
        using MPI
        using Polyester
        using HDF5
        using Aqua
        using Documenter
    end

    # macro to unlock debug mode
    DEBUG() = false

    macro DEBUG(expr, msgs)
        esc(:(if $(@__MODULE__).DEBUG() @assert($expr, $msgs...) end))
    end
    
    include("types.jl")

    include("mesh/mesh.jl")
    include("func/func.jl")

    include("misc/mpi_helpers.jl")
    include("misc/pade.jl")
    include("misc/pulay.jl")
end