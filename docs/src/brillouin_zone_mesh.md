# Brillouin zone mesh

## Constructor

```@docs
BrillouinZoneMesh
```

## Types 

```@docs
BrillouinPoint
```

```@docs
BrillouinDomain
```

```@docs
BrillouinZone
```

```@docs
LinMap
```

## Functions

```@docs
basis
```

```@docs
inv_basis
```

```@docs
bz
```

```@docs
is_inbounds(:: BrillouinPoint{N}, :: BrillouinZone{N, P}) where {N, P}
is_inbounds(:: T, :: BrillouinZone{N, P}) where {N, P, T <: AbstractVector{Float64}}
is_inbounds(:: BrillouinPoint{N}, :: Mesh{MeshPoint{BrillouinPoint{N}}, BrillouinDomain{N, P}}) where {N, P}
is_inbounds(:: T, :: Mesh{MeshPoint{BrillouinPoint{N}}, BrillouinDomain{N, P}}) where {N, P, T <: AbstractVector{Float64}}
```

```@docs
value( :: BrillouinPoint{N}) where N
```

```@docs
lin_idxs
```

```@docs
reciprocal
```

```@docs
reciprocals
```

```@docs
euclidean
```
```@docs
euclideans
```

```@docs
fold_back
```

```@docs
buffer_fold_back
```

```@docs
get_shifts
```

```@docs
MatsubaraFunctions.reverse
```

```@docs
to_Wigner_Seitz
```

```@docs
save!(:: HDF5.File, :: String, :: Mesh{MeshPoint{BrillouinPoint{N}}, BrillouinDomain{N, P}}) where{N, P}
```   

```@docs
load_mesh(:: HDF5.File, :: String, ::Val{:BrillouinZoneMesh})
```    