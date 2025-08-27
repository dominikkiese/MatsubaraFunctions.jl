# ----------------------------------------------------------------------------- #
# Symmetries Module
# ----------------------------------------------------------------------------- #

"""
	Symmetries

Defines types and functions for representing and manipulating symmetries in mesh-based functions, including the main `Symmetry` and `InitFunction` types.

Examples:
```julia
sym = Symmetry{2}(w -> (w, Operation{Float64}()))
result = sym((val1, val2))
```
"""
# load implementations
#-------------------------------------------------------------------------------#

include("operation.jl")
include("symmetry.jl")
include("symmetry_class.jl")
include("symmetry_group.jl")