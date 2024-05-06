# Index mesh

## Constructor

```@docs
IndexMesh
```

## Types
```@docs
Index
```

```@docs
IndexDomain
```


## Functions
```@docs
N(:: Mesh{MeshPoint{Index}, IndexDomain})
```


```@docs
is_inbounds(:: Index, :: Mesh{MeshPoint{Index}, IndexDomain})
is_inbounds(:: Int, :: Mesh{MeshPoint{Index}, IndexDomain})
```


```@docs
value( :: Index)
```

```@docs
values(:: Mesh{MeshPoint{Index}, IndexDomain})
```


```@docs
save!(:: HDF5.File, :: String, :: Mesh{MeshPoint{Index}, IndexDomain})
```   

```@docs
load_mesh(:: HDF5.File, :: String, ::Val{:IndexMesh})
```    