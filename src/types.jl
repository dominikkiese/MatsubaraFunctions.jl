# define particle types for dispatch
# grids can either be of fermionic or bosonic type
abstract type AbstractParticle end 

struct Fermion <: AbstractParticle end 
struct Boson   <: AbstractParticle end

# define grid types for dispatch
# grids can either be equispaced (aka linear) or coarse-grained
abstract type AbstractGrid end 

struct Linear <: AbstractGrid end 
struct Coarse <: AbstractGrid end