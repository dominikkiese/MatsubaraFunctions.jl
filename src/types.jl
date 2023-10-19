"""
    abstract type AbstractParticle
"""
abstract type AbstractParticle end 

"""
    struct Fermion <: AbstractParticle
"""
struct Fermion <: AbstractParticle end 

"""
    struct Boson <: AbstractParticle
"""
struct Boson <: AbstractParticle end

"""
    abstract type AbstractMatsubaraFrequency
"""
abstract type AbstractMatsubaraFrequency end 

"""
    abstract type AbstractMatsubaraGrid
"""
abstract type AbstractMatsubaraGrid end 

#----------------------------------------------------------------------------------------------#

export 
    AbstractParticle,
    Fermion,
    Boson,
    AbstractMatsubaraFrequency,
    AbstractMatsubaraGrid