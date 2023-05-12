# MatsubaraFrequency

A `MatsubaraFrequency` at temperature $T$ and with index $n$ can be generated using

```julia
T = 1.0
n = 5
v = MatsubaraFrequency(T, n, Fermion)
w = MatsubaraFrequency(T, n, Boson) 
```

Matsubara frequencies can be added, subtracted and their sign can be reversed, producing a new instance of `MatsubaraFrequency`.

```julia
v1 = v + v # type(v1) = :Boson
v2 = w - v # type(v2) = :Fermion
v3 = -v    # type(v3) = :Fermion 
```

# Abstract types

```@docs
AbstractParticle
```

```@docs
Fermion
```

```@docs
Boson
```

# Structs

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

```@docs
type
```