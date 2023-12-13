module MatsubaraFunctions

    using PrecompileTools

    @recompile_invalidations begin
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

    include("mesh/mesh.jl")
    include("func/func.jl")
end