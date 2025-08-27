"""
    abstract type AbstractParticle
"""

"""
    AbstractParticle

Supertype for all particle types (fermion, boson).
"""
abstract type AbstractParticle end 

"""
    struct Fermion <: AbstractParticle
"""

"""
    Fermion <: AbstractParticle

Type representing a fermionic particle.
"""
struct Fermion <: AbstractParticle end 

"""
    struct Boson <: AbstractParticle
"""

"""
    Boson <: AbstractParticle

Type representing a bosonic particle.
"""
struct Boson <: AbstractParticle end

"""
    abstract type AbstractMatsubaraFrequency
"""

"""
    AbstractMatsubaraFrequency

Supertype for Matsubara frequency types.
"""
abstract type AbstractMatsubaraFrequency end 

"""
    abstract type AbstractMatsubaraGrid
"""

"""
    AbstractMatsubaraGrid

Supertype for Matsubara grid types.
"""
abstract type AbstractMatsubaraGrid end 

#----------------------------------------------------------------------------------------------#

export 
    AbstractParticle,
    Fermion,
    Boson,
    AbstractMatsubaraFrequency,
    AbstractMatsubaraGrid