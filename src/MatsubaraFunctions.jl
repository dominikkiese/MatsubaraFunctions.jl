module MatsubaraFunctions

    using PrecompileTools

    @recompile_invalidations begin
        using Printf
        using LinearAlgebra
        using OffsetArrays
        using MPI
        using Polyester
        using HDF5
        using Documenter
        using Aqua
    end
    
    # macro to unlock debug mode
    DEBUG() = false

    macro DEBUG(expr, msgs)
        esc(:(if $(@__MODULE__).DEBUG() @assert($expr, $msgs...) end))
    end

    include("types.jl")
    include("matsubara_freq/matsubara_freq.jl")
    include("matsubara_grid/matsubara_grid.jl")
    include("matsubara_func/matsubara_func.jl")
    include("misc/mpi_helpers.jl")
    include("misc/pade.jl")
    include("misc/pulay.jl")
end