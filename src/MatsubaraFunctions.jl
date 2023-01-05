module MatsubaraFunctions

    using Test 
    using LoopVectorization
    using HDF5

    include("matsubara.jl")
    include("interpolation.jl")
    include("function.jl")
    include("io.jl")

    export mk_grid, MatsubaraFunction, add!, subtract!, mult!, save_matsubara_function!, load_matsubara_function
end
