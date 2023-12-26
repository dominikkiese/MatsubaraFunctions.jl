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



# abstract types for meshes
#-------------------------------------------------------------------------------#

"""
    abstract type AbstractMeshPoint

AbstractMeshPoint type
"""
abstract type AbstractMeshPoint end

"""
    abstract type AbstractValue

AbstractValue type
"""
abstract type AbstractValue end

"""
    abstract type AbstractMesh

AbstractMesh type
"""
abstract type AbstractMesh end



#----------------------------------------------------------------------------------------------#

export 
    AbstractParticle,
    Fermion,
    Boson,
    AbstractMeshPoint,
    AbstractValue,
    AbstractMesh