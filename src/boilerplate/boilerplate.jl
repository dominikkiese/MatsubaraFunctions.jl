@generated function _mesh_indices(f :: MeshFunction{DD, Q, MT, AT}, x :: Vararg{Union{MeshPoint, <: AbstractValue, Int, UnitRange, Colon, Float64, <: AbstractVector{Float64}}, DD}
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}
    
    return Meta.parse("("*prod(["mesh_index(x[$i], meshes(f, $i))," for i in 1 : DD])*")")
end

@generated function _mesh_indices_bc(f :: MeshFunction{DD, Q, MT, AT}, x :: Vararg{Union{MeshPoint, <: AbstractValue, Int, UnitRange, Colon, Float64, <: AbstractVector{Float64}}, DD}
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    return Meta.parse("("*prod(["mesh_index_bc(x[$i], meshes(f, $i))," for i in 1 : DD])*")")
end

@generated function _all_inbounds(f :: MeshFunction{DD, Q, MT, AT}, p :: Vararg{Union{MeshPoint, <: AbstractValue, Int, UnitRange, Colon, Float64, <: AbstractVector{Float64}}, DD}
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    return Meta.parse("all(("*prod(["is_inbounds(p[$i], meshes(f, $i))," for i in 1 : DD])*"))")
end

@generated function _all_inbounds_bc(f :: MeshFunction{DD, Q, MT, AT}, p :: Vararg{Union{MeshPoint, <: AbstractValue, Int, UnitRange, Colon, Float64, <: AbstractVector{Float64}}, DD}
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    return Meta.parse("all(("*prod(["is_inbounds_bc(p[$i], meshes(f, $i))," for i in 1 : DD])*"))")
end

@generated function _get_params(f :: MeshFunction{DD, Q, MT, AT}, p :: Vararg{Union{MeshPoint, <: AbstractValue, Int, Float64, <: AbstractVector{Float64}}, DD}
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    return Meta.parse("("*prod(["InterpolationParam(p[$i], meshes(f, $i))," for i in 1 : DD])*")")
end

"""
    function to_meshes(f :: MeshFunction{DD, Q, MT, AT}, cidx :: CartesianIndex{DD}
        ) :: NTuple{DD, Union{MeshPoint}} where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

Returns mesh points
"""
@generated function to_meshes(f :: MeshFunction{DD, Q, MT, AT}, cidx :: CartesianIndex{DD}
    ) :: NTuple{DD, Union{MeshPoint}} where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    return Meta.parse("("*prod(["meshes(f, $i)[cidx[$i]]," for i in 1 : DD])*")")
end