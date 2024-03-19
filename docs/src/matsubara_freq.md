# MatsubaraFrequency

A `MatsubaraFrequency{PT}` of particle type `PT` at temperature $T$ and with Matsubara index $n$ can be generated using

```julia
T = 1.0
n = 5
v = MatsubaraFrequency(T, n, Fermion) # v = \pi T (2 n + 1) 
w = MatsubaraFrequency(T, n, Boson)   # w = \pi T (2 n    )
```

Matsubara frequencies can be added, subtracted and their sign can be reversed, producing a new instance of `MatsubaraFrequency`, 
potentially of different particle type.

```julia
v1 = v + v # typeof(v1) = MatsubaraFrequency{Boson}
v2 = w - v # typeof(v2) = MatsubaraFrequency{Fermion}
v3 = -v    # typeof(v3) = MatsubaraFrequency{Fermion}
```

# Types

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
AbstractMatsubaraFrequency
```

```@docs
MatsubaraFrequency
```

# Functions

```@docs
temperature
```

```@docs
value
```

```@docs
index
```