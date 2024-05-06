# Matsubara mesh

## Constructor

```@docs
MatsubaraMesh
```

## Types

```@docs
AbstractParticle
```

```@docs
Fermion
```

```@docs
Boson
```

```@docs
MatsubaraFrequency
```

```@docs
MatsubaraDomain
```

## Functions

```@docs
temperature
```

```@docs
N(:: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) where {PT <: AbstractParticle}
```

```@docs
index(:: MatsubaraFrequency{PT}) where {PT <: AbstractParticle}
```

```@docs
indices(:: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) where PT<:AbstractParticle
```

```@docs
first_index
```

```@docs
last_index
```

```@docs
is_inbounds(:: MatsubaraFrequency{PT}, :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) where {PT <: AbstractParticle}
is_inbounds(:: Float64, :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) where {PT <: AbstractParticle}
```

```@docs
value(:: MatsubaraFrequency{PT}) where {PT <: AbstractParticle}
```

```@docs
values(:: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) where {PT <: AbstractParticle}
```

```@docs
first_value
```

```@docs
last_value
```

```@docs
save!(:: HDF5.File, :: String, :: Mesh{MeshPoint{MatsubaraFrequency{Fermion}}, MatsubaraDomain})
save!(:: HDF5.File, :: String, :: Mesh{MeshPoint{MatsubaraFrequency{Boson}}, MatsubaraDomain})
```   

```@docs
load_mesh(:: HDF5.File, :: String, ::Val{:MatsubaraMesh})
```    