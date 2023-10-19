module MatsubaraFunctions

    using PrecompileTools

    @recompile_invalidations begin
        using Printf
        using LinearAlgebra
        using MPI
        using Polyester
        using HDF5
        using Documenter
        using Aqua
    end

    # macro for sanity checks
    sanity_checks() = true 

    macro check(expr, msgs)
        esc(:(if $(@__MODULE__).sanity_checks() @assert($expr, $msgs...) end))
    end

    include("types.jl")
    include("matsubara_freq/matsubara_freq.jl")
    include("matsubara_grid/matsubara_grid.jl")
    include("matsubara_func/matsubara_func.jl")
    include("misc/mpi_helpers.jl")
    include("misc/pade.jl")
    include("misc/pulay.jl")
end