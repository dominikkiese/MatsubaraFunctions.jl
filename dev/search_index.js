var documenterSearchIndex = {"docs":
[{"location":"#MatsubaraFunctions.jl-Documentation","page":"Home","title":"MatsubaraFunctions.jl Documentation","text":"","category":"section"},{"location":"#Abstract-Types","page":"Home","title":"Abstract Types","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"AbstractParticle","category":"page"},{"location":"#MatsubaraFunctions.AbstractParticle","page":"Home","title":"MatsubaraFunctions.AbstractParticle","text":"abstract type AbstractParticle\n\nAbstractParticle type\n\n\n\n\n\n","category":"type"},{"location":"","page":"Home","title":"Home","text":"Fermion","category":"page"},{"location":"#MatsubaraFunctions.Fermion","page":"Home","title":"MatsubaraFunctions.Fermion","text":"struct Fermion <: AbstractParticle\n\nFermionic particle type, used for MatsubaraFrequency and MatsubaraGrid constructors\n\n\n\n\n\n","category":"type"},{"location":"","page":"Home","title":"Home","text":"Boson","category":"page"},{"location":"#MatsubaraFunctions.Boson","page":"Home","title":"MatsubaraFunctions.Boson","text":"struct Boson <: AbstractParticle\n\nBosonic particle type, used for MatsubaraFrequency and MatsubaraGrid constructors\n\n\n\n\n\n","category":"type"},{"location":"#Structs","page":"Home","title":"Structs","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"MatsubaraFrequency","category":"page"},{"location":"#MatsubaraFunctions.MatsubaraFrequency","page":"Home","title":"MatsubaraFunctions.MatsubaraFrequency","text":"struct MatsubaraFrequency\n\nMatsubaraFrequency type with fields:\n\nT    :: Float64 : physical temperature\nval  :: Float64 : position on the imaginary axis\nidx  :: Int64   : Matsubara index\ntype :: Symbol  : particle type\n\nExamples:\n\n# construction\nT   = 1.0\nidx = 5\nv   = MatsubaraFrequency(T, idx, Fermion)\nw   = MatsubaraFrequency(T, idx, Boson) \n\n# usage\nw1 = v + v # type(v1) = :Boson\nv2 = w - v # type(v2) = :Fermion\nv3 = -v    # type(v3) = :Fermion \n\n\n\n\n\n","category":"type"},{"location":"#Functions","page":"Home","title":"Functions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"temperature","category":"page"},{"location":"#MatsubaraFunctions.temperature","page":"Home","title":"MatsubaraFunctions.temperature","text":"function temperature(\n    w :: MatsubaraFrequency\n    ) :: Float64\n\nReturns w.T\n\n\n\n\n\nfunction temperature(\n    grid :: MatsubaraGrid\n    )    :: Float64\n\nReturns grid.T\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"value","category":"page"},{"location":"#MatsubaraFunctions.value","page":"Home","title":"MatsubaraFunctions.value","text":"function value(\n    w :: MatsubaraFrequency\n    ) :: Float64\n\nReturns w.val\n\n\n\n\n\n","category":"function"}]
}
