module MatsubaraFunctions

    using Test 
    using LoopVectorization

    include("matsubara.jl")
    include("interpolation.jl")
    include("function.jl")

    export mk_grid, MatsubaraFunction, add!, subtract!, mult!
end
