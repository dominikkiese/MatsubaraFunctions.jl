"""
    function integrate(f :: Function, ms :: Vararg{<: AbstractMesh})

Integrate the function `f` over the meshes `ms`
"""
function integrate(f :: Function, ms :: Vararg{<: AbstractMesh})
    vol = prod(volume_element(m) for m in ms)
    s   = ThreadsX.sum(f(iter...) for iter in collect(Iterators.product(ms...)))
    return vol * s
end

"""
    function integrate(f :: MeshFunction{DD, Q, MT, AT}) :: Q where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

Integrate the mesh function `f` over the mesh defined by the meshes in `MT`
"""
function integrate(f :: MeshFunction{DD, Q, MT, AT}) :: Q where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}
    return integrate((x...) -> f(x...), meshes(f)...)
end

export integrate