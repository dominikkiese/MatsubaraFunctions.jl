var documenterSearchIndex = {"docs":
[{"location":"matsubara_freq/#MatsubaraFrequency","page":"MatsubaraFrequency","title":"MatsubaraFrequency","text":"","category":"section"},{"location":"matsubara_freq/","page":"MatsubaraFrequency","title":"MatsubaraFrequency","text":"A MatsubaraFrequency at temperature T and with index n can be generated using","category":"page"},{"location":"matsubara_freq/","page":"MatsubaraFrequency","title":"MatsubaraFrequency","text":"T = 1.0\nn = 5\nv = MatsubaraFrequency(T, n, Fermion)\nw = MatsubaraFrequency(T, n, Boson) ","category":"page"},{"location":"matsubara_freq/","page":"MatsubaraFrequency","title":"MatsubaraFrequency","text":"Matsubara frequencies can be added, subtracted and their sign can be reversed, producing a new instance of MatsubaraFrequency.","category":"page"},{"location":"matsubara_freq/","page":"MatsubaraFrequency","title":"MatsubaraFrequency","text":"v1 = v + v # type(v1) = :Boson\nv2 = w - v # type(v2) = :Fermion\nv3 = -v    # type(v3) = :Fermion ","category":"page"},{"location":"matsubara_freq/#Abstract-types","page":"MatsubaraFrequency","title":"Abstract types","text":"","category":"section"},{"location":"matsubara_freq/","page":"MatsubaraFrequency","title":"MatsubaraFrequency","text":"AbstractParticle","category":"page"},{"location":"matsubara_freq/#MatsubaraFunctions.AbstractParticle","page":"MatsubaraFrequency","title":"MatsubaraFunctions.AbstractParticle","text":"abstract type AbstractParticle\n\nAbstractParticle type\n\n\n\n\n\n","category":"type"},{"location":"matsubara_freq/","page":"MatsubaraFrequency","title":"MatsubaraFrequency","text":"Fermion","category":"page"},{"location":"matsubara_freq/#MatsubaraFunctions.Fermion","page":"MatsubaraFrequency","title":"MatsubaraFunctions.Fermion","text":"struct Fermion <: AbstractParticle\n\nFermionic particle type, used for MatsubaraFrequency and MatsubaraGrid constructors\n\n\n\n\n\n","category":"type"},{"location":"matsubara_freq/","page":"MatsubaraFrequency","title":"MatsubaraFrequency","text":"Boson","category":"page"},{"location":"matsubara_freq/#MatsubaraFunctions.Boson","page":"MatsubaraFrequency","title":"MatsubaraFunctions.Boson","text":"struct Boson <: AbstractParticle\n\nBosonic particle type, used for MatsubaraFrequency and MatsubaraGrid constructors\n\n\n\n\n\n","category":"type"},{"location":"matsubara_freq/#Structs","page":"MatsubaraFrequency","title":"Structs","text":"","category":"section"},{"location":"matsubara_freq/","page":"MatsubaraFrequency","title":"MatsubaraFrequency","text":"MatsubaraFrequency","category":"page"},{"location":"matsubara_freq/#MatsubaraFunctions.MatsubaraFrequency","page":"MatsubaraFrequency","title":"MatsubaraFunctions.MatsubaraFrequency","text":"struct MatsubaraFrequency\n\nMatsubaraFrequency type with fields:\n\nT    :: Float64 : physical temperature\nval  :: Float64 : position on the imaginary axis\nidx  :: Int64   : Matsubara index\ntype :: Symbol  : particle type\n\nExamples:\n\n# construction\nT   = 1.0\nidx = 5\nv   = MatsubaraFrequency(T, idx, Fermion)\nw   = MatsubaraFrequency(T, idx, Boson) \n\n# usage\nw1 = v + v # type(v1) = :Boson\nv2 = w - v # type(v2) = :Fermion\nv3 = -v    # type(v3) = :Fermion \n\n\n\n\n\n","category":"type"},{"location":"matsubara_freq/#Functions","page":"MatsubaraFrequency","title":"Functions","text":"","category":"section"},{"location":"matsubara_freq/","page":"MatsubaraFrequency","title":"MatsubaraFrequency","text":"temperature","category":"page"},{"location":"matsubara_freq/#MatsubaraFunctions.temperature","page":"MatsubaraFrequency","title":"MatsubaraFunctions.temperature","text":"function temperature(\n    w :: MatsubaraFrequency\n    ) :: Float64\n\nReturns w.T\n\n\n\n\n\nfunction temperature(\n    grid :: MatsubaraGrid\n    )    :: Float64\n\nReturns grid.T\n\n\n\n\n\n","category":"function"},{"location":"matsubara_freq/","page":"MatsubaraFrequency","title":"MatsubaraFrequency","text":"value","category":"page"},{"location":"matsubara_freq/#MatsubaraFunctions.value","page":"MatsubaraFrequency","title":"MatsubaraFunctions.value","text":"function value(\n    w :: MatsubaraFrequency\n    ) :: Float64\n\nReturns w.val\n\n\n\n\n\n","category":"function"},{"location":"matsubara_freq/","page":"MatsubaraFrequency","title":"MatsubaraFrequency","text":"index","category":"page"},{"location":"matsubara_freq/#MatsubaraFunctions.index","page":"MatsubaraFrequency","title":"MatsubaraFunctions.index","text":"function index(\n    w :: MatsubaraFrequency\n    ) :: Int64\n\nReturns w.idx\n\n\n\n\n\n","category":"function"},{"location":"matsubara_freq/","page":"MatsubaraFrequency","title":"MatsubaraFrequency","text":"type","category":"page"},{"location":"matsubara_freq/#MatsubaraFunctions.type","page":"MatsubaraFrequency","title":"MatsubaraFunctions.type","text":"function type(\n    w :: MatsubaraFrequency\n    ) :: Symbol\n\nReturns w.type\n\n\n\n\n\nfunction type(\n    grid :: MatsubaraGrid\n    )    :: Symbol\n\nReturns grid.type\n\n\n\n\n\n","category":"function"},{"location":"matsubara_grid/#MatsubaraGrid","page":"MatsubaraGrid","title":"MatsubaraGrid","text":"","category":"section"},{"location":"matsubara_grid/","page":"MatsubaraGrid","title":"MatsubaraGrid","text":"A MatsubaraGrid is a sorted (symmetric) set of MatsubaraFrequency objects and can be constructed by","category":"page"},{"location":"matsubara_grid/","page":"MatsubaraGrid","title":"MatsubaraGrid","text":"T  = 1.0\nN  = 128\ng1 = MatsubaraGrid(T, N, Fermion) # no. frequencies is 2N\ng2 = MatsubaraGrid(T, N, Boson)   # no. frequencies is 2N - 1","category":"page"},{"location":"matsubara_grid/","page":"MatsubaraGrid","title":"MatsubaraGrid","text":"MatsubaraGrid instances are iterable","category":"page"},{"location":"matsubara_grid/","page":"MatsubaraGrid","title":"MatsubaraGrid","text":"T = 1.0\nN = 128\ng = MatsubaraGrid(T, N, Fermion)\n\nfor v in g\n  println(value(v)) \n  println(index(v))\nend","category":"page"},{"location":"matsubara_grid/","page":"MatsubaraGrid","title":"MatsubaraGrid","text":"and can be evaluated using either a MatsubaraFrequency or Float64. As long as the input argument is in bounds, this will return the corresponding linear index of the grid in the former case and the linear index of the closest frequency in the latter case ","category":"page"},{"location":"matsubara_grid/","page":"MatsubaraGrid","title":"MatsubaraGrid","text":"T   = 1.0\nN   = 128\ng   = MatsubaraGrid(T, N, Fermion)\nidx = rand(1 : length(g))\n@assert g(g[idx]) == idx \n@assert g(value(g[idx])) == idx ","category":"page"},{"location":"matsubara_grid/","page":"MatsubaraGrid","title":"MatsubaraGrid","text":"MatsubaraGrid objects can be saved in HDF5 file format as","category":"page"},{"location":"matsubara_grid/","page":"MatsubaraGrid","title":"MatsubaraGrid","text":"``julia using MatsubaraFunctions  using HDF5","category":"page"},{"location":"matsubara_grid/","page":"MatsubaraGrid","title":"MatsubaraGrid","text":"file = h5open(\"test.h5\", \"w\") T    = 1.0 N    = 128 g    = MatsubaraGrid(T, N, Fermion)","category":"page"},{"location":"matsubara_grid/","page":"MatsubaraGrid","title":"MatsubaraGrid","text":"savematsubaragrid!(file, \"grid\", g)  gp = loadmatsubaragrid(file, \"grid\") close(file)","category":"page"},{"location":"matsubara_grid/","page":"MatsubaraGrid","title":"MatsubaraGrid","text":"\n# Structs\n","category":"page"},{"location":"matsubara_grid/","page":"MatsubaraGrid","title":"MatsubaraGrid","text":"@docs MatsubaraGrid","category":"page"},{"location":"matsubara_grid/","page":"MatsubaraGrid","title":"MatsubaraGrid","text":"\n# Functions\n","category":"page"},{"location":"matsubara_grid/","page":"MatsubaraGrid","title":"MatsubaraGrid","text":"@docs index_range","category":"page"},{"location":"matsubara_grid/","page":"MatsubaraGrid","title":"MatsubaraGrid","text":"","category":"page"},{"location":"matsubara_grid/","page":"MatsubaraGrid","title":"MatsubaraGrid","text":"@docs is_inbounds","category":"page"},{"location":"matsubara_grid/","page":"MatsubaraGrid","title":"MatsubaraGrid","text":"","category":"page"},{"location":"matsubara_grid/","page":"MatsubaraGrid","title":"MatsubaraGrid","text":"@docs info","category":"page"},{"location":"matsubara_grid/","page":"MatsubaraGrid","title":"MatsubaraGrid","text":"","category":"page"},{"location":"matsubara_grid/","page":"MatsubaraGrid","title":"MatsubaraGrid","text":"@docs savematsubaragrid!","category":"page"},{"location":"matsubara_grid/","page":"MatsubaraGrid","title":"MatsubaraGrid","text":"","category":"page"},{"location":"matsubara_grid/","page":"MatsubaraGrid","title":"MatsubaraGrid","text":"@docs loadmatsubaragrid ```    ","category":"page"},{"location":"matsubara_func/#MatsubaraFunction","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"","category":"section"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"A MatsubaraFunction is a collection of MatsubaraGrid instances together with an associated tensor structure G_i_1i_n for each point (omega_1  omega_m) in the cartesian product of the grids.  Some possible constructors are","category":"page"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"T = 1.0\nN = 128\ng = MatsubaraGrid(T, N, Fermion)\n\n# 1D grid, rank 1 tensor with index dimension 1 (scalar valued)\nf1_complex = MatsubaraFunction(g, 1) \nf1_real    = MatsubaraFunction(g, 1, Float64) \ninfo(f1_complex)\ninfo(f1_real)\n\n# 1D grid, rank 1 tensor with index dimension 5 (vector valued)\nf2_complex = MatsubaraFunction(g, 5) \nf2_real    = MatsubaraFunction(g, 5, Float64) \ninfo(f2_complex)\ninfo(f2_real)\n\n# 1D grid, rank 2 tensor with index dimension 5 (matrix valued)\nf3_complex = MatsubaraFunction(g, (5, 5)) \nf3_real    = MatsubaraFunction(g, (5, 5), Float64) \ninfo(f3_complex)\ninfo(f3_real)\n\n# 2D grid, rank 2 tensor with index dimension 5 (matrix valued)\nf4_complex = MatsubaraFunction((g, g), (5, 5)) \nf4_real    = MatsubaraFunction((g, g), (5, 5), Float64) \ninfo(f4_complex)\ninfo(f4_real)","category":"page"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"There are two possible ways to access the data of a MatsubaraFunction, using either the bracket [] or the parenthesis () operator. The former can be used together with a set of linear indices or with a combination of MatsubaraFrequency objects and linear indices (for the tensor structure). It will return the value of the function precisely for the given arguments. () allows to substitute Float64 for the frequency arguments, in which case a multilinear interpolation is performed. In addition, () is well defined even for out of bounds access, using either a custom boundary condition or, for 1D grids, polynomial extrapolation.","category":"page"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"ξ = 0.5\nT = 1.0\nN = 128\ng = MatsubaraGrid(T, N, Fermion)\nf = MatsubaraFunction(g, 1)\n\nfor v in g\n    f[v, 1] = 1.0 / (im * value(v) - ξ)\nend \n\n# access MatsubaraFunction data\nprintln(f[v, 1])        # fast data access, throws error if v is out of bounds\nprintln(f(v, 1))        # fast data access, defined even if v is out of bounds\nprintln(f(value(v), 1)) # slow data access, uses interpolation \n\n# fallback methods for out of bounds access\nvp = MatsubaraFrequency(T, 256, Fermion)\nprintln(f(vp, 1))                                  # default x -> 0.0\nprintln(f(vp, 1; bc = x -> 1.0))                   # custom boundary condition x -> 1.0\nprintln(f(vp, 1; bc = x -> 1.0 / im / value(x)))   # custom boundary condition x -> 1.0 / im / value(x)\nprintln(f(value(vp), 1; bc = x -> 1.0 / im / x))   # custom boundary condition x -> 1.0 / im / x\nprintln(f(vp, 1; extrp = (true, ComplexF64(0.0)))) # polynomial extrapolation in 1D, constant term set to 0.0","category":"page"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"MatsubaraFunction objects can be saved in HDF5 file format as","category":"page"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"using MatsubaraFunctions \nusing HDF5\n\nfile = h5open(\"test.h5\", \"w\")\nT    = 1.0\nN    = 128\ng    = MatsubaraGrid(T, N, Fermion)\nf    = MatsubaraFunction(g, 1)\n\nsave_matsubara_function!(file, \"func\", f)\nfp = load_matsubara_function(file, \"func\")\nclose(file)","category":"page"},{"location":"matsubara_func/#Advanced-Usage:-Matsubara-Sums","page":"MatsubaraFunction","title":"Advanced Usage: Matsubara Sums","text":"","category":"section"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"For MatsubaraFunction objects G_i_1  i_n(iomega) defined on 1D grids, we export the function sum_me, which computes the series Sigma_m G_i_1  i_n(iomega^m) e^iomega^m 0^+ for m in mathbbZ using tail fits of G together with analytic formulas for summations of the form Sigma_m frac1(iomega^m)^alphae^iomega^m 0^+ with alpha in mathbbN. This, however, requires G to be representable by a Laurent series in an elongated annulus about the imaginary axis.","category":"page"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"ξ = 0.5\nT = 1.0\nN = 128\ng = MatsubaraGrid(T, N, Fermion)\nf = MatsubaraFunction(g, 1)\n\nfor v in g\n    f[v, 1] = 1.0 / (im * value(v) - ξ)\nend \n\n# evaluate the series and compare to analytic result\nρ(x, T) = 1.0 / (exp(x / T) + 1.0)\nprintln(abs(sum_me(f, ComplexF64(0.0), 1) - (ρ(+ξ, T) - 1.0)))","category":"page"},{"location":"matsubara_func/#Advanced-Usage:-Automated-Symmetry-Reduction","page":"MatsubaraFunction","title":"Advanced Usage: Automated Symmetry Reduction","text":"","category":"section"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"In many cases, the numerical effort of computing functions in the Matsubara domain can be drastically reduced by the use of symmetries. For one-particle Green's functions G_i_1 i_2(iomega), for example, hermicity of the Hamiltonian dictates that G_i_1 i_2(iomega) = G^star_i_2 i_1(-iomega), relating positive and negative Matsubara frequencies. We offer an automated way to compute the set of irreducible (i.e. unrelatable by symmetries) MatsubaraFunction components, as is illustrated in the following example","category":"page"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"ξ = 0.5\nT = 1.0\nN = 128\ng = MatsubaraGrid(T, N, Fermion)\nf = MatsubaraFunction(g, 1)\n\nfor v in g\n    f[v, 1] = 1.0 / (im * value(v) - ξ)\nend \n\n# complex conjugation acting on Green's function\nfunction conj(\n    w :: Tuple{MatsubaraFrequency},\n    x :: Tuple{Int64}\n    ) :: Tuple{Tuple{MatsubaraFrequency}, Tuple{Int64}, MatsubaraOperation}\n\n    return (-w[1],), (x[1],), MatsubaraOperation(false, true)\nend \n\n# compute the symmetry group \nSG = MatsubaraSymmetryGroup([MatsubaraSymmetry{1, 1}(conj)], f)\n\n# obtain another Green's function by symmetrization\nfunction init(\n    w :: Tuple{MatsubaraFrequency},\n    x :: Tuple{Int64}\n    ) :: ComplexF64\n\n    return f[w, x...]\nend \n\nInitFunc = MatsubaraInitFunction{1, 1, ComplexF64}(init)\nh        = MatsubaraFunction(g, 1)\nSG(h, InitFunc)\n@assert h == f","category":"page"},{"location":"matsubara_func/#Advanced-Usage:-MPI-Helpers","page":"MatsubaraFunction","title":"Advanced Usage: MPI Helpers","text":"","category":"section"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"To simplify the parallelization of algorithms involving MatsubaraFunction instances, we export some useful methods based on the MPI.jl wrapper. For further information on how to set up MPI with Julia see https://github.com/JuliaParallel/MPI.jl.","category":"page"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"using MatsubaraFunctions \nusing MPI \n\nMPI.Init()\nmpi_info()\nmpi_println(\"I print on main.\")\nismain = mpi_ismain() # ismain = true if rank is 0\n\nT = 1.0\nN = 128\ng = MatsubaraGrid(T, N, Fermion)\nf = MatsubaraFunction(g, 1)\n\n# simple loop parallelization for UnitRange\nfor vidx in mpi_split(1 : length(g))\n  comm = MPI.COMM_WORLD \n  println(\"My rank is $(MPI.Comm_rank(comm)): $(vidx)\")\nend\n\n# simple (+) allreduce\nmpi_allreduce!(f)","category":"page"},{"location":"matsubara_func/#Structs","page":"MatsubaraFunction","title":"Structs","text":"","category":"section"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"MatsubaraFunction","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.MatsubaraFunction","page":"MatsubaraFunction","title":"MatsubaraFunctions.MatsubaraFunction","text":"struct MatsubaraFunction{GD, SD, DD, Q <: Number}\n\nMatsubaraFunction type with fields:\n\ngrids :: NTuple{GD, MatsubaraGrid} : collection of MatsubaraGrid\nshape :: NTuple{SD, Int64}         : shape of the tensor structure on every grid point\ndata  :: Array{Q, DD}              : data array\n\nExamples:\n\n# construction\nT = 1.0\nN = 128\ng = MatsubaraGrid(T, N, Fermion)\n\n# 1D grid, rank 1 tensor with index dimension 1 (scalar valued)\nf1_complex = MatsubaraFunction(g, 1)                    # complex valued (default)\nf1_real    = MatsubaraFunction(g, 1, Float64)           # other data type\n\n# 1D grid, rank 1 tensor with index dimension 5 (vector valued)\nf2_complex = MatsubaraFunction(g, 5)                    # complex valued (default)\nf2_real    = MatsubaraFunction(g, 5, Float64)           # other data type\n\n# 1D grid, rank 2 tensor with index dimension 5 (matrix valued)\nf3_complex = MatsubaraFunction(g, (5, 5))               # complex valued (default)\nf3_real    = MatsubaraFunction(g, (5, 5), Float64)      # other data type\n\n# 2D grid, rank 2 tensor with index dimension 5 (matrix valued)\nf4_complex = MatsubaraFunction((g, g), (5, 5))          # complex valued (default)\nf4_real    = MatsubaraFunction((g, g), (5, 5), Float64) # other data type\n\n# usage \nξ = 0.5\nf = MatsubaraFunction(g, 1)\n\nfor v in g\n    f[v, 1] = 1.0 / (im * value(v) - ξ)\nend \n\n# access MatsubaraFunction data\nprintln(f[v, 1])        # fast data access, throws error if v is out of bounds\nprintln(f(v, 1))        # fast data access, defined even if v is out of bounds\nprintln(f(value(v), 1)) # slow data access, uses interpolation \n\n# fallback methods for out of bounds access\nvp = MatsubaraFrequency(T, 256, Fermion)\nprintln(f(vp, 1))                                  # default x -> 0.0\nprintln(f(vp, 1; bc = x -> 1.0))                   # custom boundary condition x -> 1.0\nprintln(f(vp, 1; bc = x -> 1.0 / im / value(x)))   # custom boundary condition x -> 1.0 / im / value(x)\nprintln(f(value(vp), 1; bc = x -> 1.0 / im / x))   # custom boundary condition x -> 1.0 / im / x\nprintln(f(vp, 1; extrp = (true, ComplexF64(0.0)))) # polynomial extrapolation in 1D, constant term set to 0.0\n\n\n\n\n\n","category":"type"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"MatsubaraOperation","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.MatsubaraOperation","page":"MatsubaraFunction","title":"MatsubaraFunctions.MatsubaraOperation","text":"struct MatsubaraOperation\n\nMatsubaraOperation type with fields:\n\nsgn :: Bool : change sign?\ncon :: Bool : do complex conjugation?\n\n\n\n\n\n","category":"type"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"MatsubaraSymmetry","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.MatsubaraSymmetry","page":"MatsubaraFunction","title":"MatsubaraFunctions.MatsubaraSymmetry","text":"struct MatsubaraSymmetry{GD, SD}\n\nMatsubaraSymmetry type with fields:\n\nf :: FunctionWrappers.FunctionWrapper{Tuple{NTuple{GD, MatsubaraFrequency}, NTuple{SD, Int64}, MatsubaraOperation}, Tuple{NTuple{GD, MatsubaraFrequency}, NTuple{SD, Int64}}} :\n\nMatsubaraSymmetry takes grid coordinates and tensor indices as input and returns a new set of coordinates and indices together with a MatsubaraOperation\n\n\n\n\n\n","category":"type"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"MatsubaraSymmetryGroup","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.MatsubaraSymmetryGroup","page":"MatsubaraFunction","title":"MatsubaraFunctions.MatsubaraSymmetryGroup","text":"MatsubaraSymmetryGroup\n\nMatsubaraSymmetryGroup type with fields:\n\nclasses :: Vector{Vector{Tuple{Int64, MatsubaraOperation}}} : list of collections of symmetry equivalent elements\n\nExamples:\n\n# a simple Green's function\nξ = 0.5\nT = 1.0\nN = 128\ng = MatsubaraGrid(T, N, Fermion)\nf = MatsubaraFunction(g, 1)\n\nfor v in g\n    f[v, 1] = 1.0 / (im * value(v) - ξ)\nend \n\n# complex conjugation acting on Green's function\nfunction conj(\n    w :: Tuple{MatsubaraFrequency},\n    x :: Tuple{Int64}\n    ) :: Tuple{Tuple{MatsubaraFrequency}, Tuple{Int64}, MatsubaraOperation}\n\n    return (-w[1],), (x[1],), MatsubaraOperation(false, true)\nend \n\n# compute the symmetry group \nSG = MatsubaraSymmetryGroup([MatsubaraSymmetry{1, 1}(conj)], f)\n\n# obtain another Green's function by symmetrization\nfunction init(\n    w :: Tuple{MatsubaraFrequency},\n    x :: Tuple{Int64}\n    ) :: ComplexF64\n\n    return f[w, x...]\nend \n\nInitFunc = MatsubaraInitFunction{1, 1, ComplexF64}(init)\nh        = MatsubaraFunction(g, 1)\nSG(h, InitFunc)\n@assert h == f\n\n\n\n\n\n","category":"type"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"MatsubaraInitFunction","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.MatsubaraInitFunction","page":"MatsubaraFunction","title":"MatsubaraFunctions.MatsubaraInitFunction","text":"struct MatsubaraInitFunction{GD, SD, Q <: Number}\n\nMatsubaraInitFunction type with fields:\n\nf :: FunctionWrappers.FunctionWrapper{Q, Tuple{NTuple{GD, MatsubaraFrequency}, NTuple{SD, Int64}}} \n\nMatsubaraInitFunction takes grid coordinates and tensor indices as input and returns value of type Q\n\n\n\n\n\n","category":"type"},{"location":"matsubara_func/#Functions","page":"MatsubaraFunction","title":"Functions","text":"","category":"section"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"grids_shape","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.grids_shape","page":"MatsubaraFunction","title":"MatsubaraFunctions.grids_shape","text":"function grids_shape(\n    f :: MatsubaraFunction{GD, SD, DD, Q}\n    ) :: NTuple{GD, Int64} where {GD, SD, DD, Q <: Number}\n\nReturns length of grids\n\n\n\n\n\nfunction grids_shape(\n    f   :: MatsubaraFunction{GD, SD, DD, Q},\n    idx :: Int64\n    )   :: Int64 where {GD, SD, DD, Q <: Number}\n\nReturns length of f.grids[idx]\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"shape","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.shape","page":"MatsubaraFunction","title":"MatsubaraFunctions.shape","text":"function shape(\n    f :: MatsubaraFunction{GD, SD, DD, Q}\n    ) :: NTuple{SD, Int64} where {GD, SD, DD, Q <: Number}\n\nReturns f.shape\n\n\n\n\n\nfunction shape(\n    f   :: MatsubaraFunction{GD, SD, DD, Q},\n    idx :: Int64\n    )   :: Int64 where {GD, SD, DD, Q <: Number}\n\nReturns f.shape[idx]\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"data_shape","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.data_shape","page":"MatsubaraFunction","title":"MatsubaraFunctions.data_shape","text":"function data_shape(\n    f :: MatsubaraFunction{GD, SD, DD, Q}\n    ) :: NTuple{DD, Int64} where {GD, SD, DD, Q <: Number}\n\nReturns shape of f.data\n\n\n\n\n\nfunction data_shape(\n    f   :: MatsubaraFunction{GD, SD, DD, Q},\n    idx :: Int64\n    )   :: Int64 where {GD, SD, DD, Q <: Number}\n\nReturns length of dimension idx of f.data\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"absmax","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.absmax","page":"MatsubaraFunction","title":"MatsubaraFunctions.absmax","text":"function absmax(\n    f :: MatsubaraFunction{GD, SD, DD, Q}\n    ) :: Float64 where {GD, SD, DD, Q <: Number}\n\nReturns largest element of f.data (in absolute terms)\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"argmax","category":"page"},{"location":"matsubara_func/#Base.argmax","page":"MatsubaraFunction","title":"Base.argmax","text":"function Base.:argmax(\n    f :: MatsubaraFunction{GD, SD, DD, Q}\n    ) :: CartesianIndex{DD} where {GD, SD, DD, Q <: Number}\n\nReturns position of largest element of f.data (in absolute terms)\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"mpi_split","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.mpi_split","page":"MatsubaraFunction","title":"MatsubaraFunctions.mpi_split","text":"function mpi_split(\n    r :: UnitRange{Int64}\n    ) :: UnitRange{Int64}\n\nSplits UnitRange evenly among available MPI ranks (including main).  Can, for example, be used to parallelize loops.\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"mpi_allreduce!","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.mpi_allreduce!","page":"MatsubaraFunction","title":"MatsubaraFunctions.mpi_allreduce!","text":"function mpi_allreduce!(\n    f :: MatsubaraFunction{GD, SD, DD, Q}\n    ) :: Nothing where {GD, SD, DD, Q <: Number}\n\nInplace MPI reduction (+) for MatsubaraFunction\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"mpi_ismain","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.mpi_ismain","page":"MatsubaraFunction","title":"MatsubaraFunctions.mpi_ismain","text":"mpi_ismain() :: Bool\n\nReturns true for MPI rank 0\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"mpi_println","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.mpi_println","page":"MatsubaraFunction","title":"MatsubaraFunctions.mpi_println","text":"function mpi_println(\n    s :: String\n    ) :: Nothing\n\nPrint string s on MPI rank 0\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"mpi_info","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.mpi_info","page":"MatsubaraFunction","title":"MatsubaraFunctions.mpi_info","text":"mpi_info() :: Nothing\n\nPrint information about available resources\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"mpi_barrier","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.mpi_barrier","page":"MatsubaraFunction","title":"MatsubaraFunctions.mpi_barrier","text":"mpi_barrier() :: Nothing\n\nPlace synchronization barrier for MPI ranks\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"add","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.add","page":"MatsubaraFunction","title":"MatsubaraFunctions.add","text":"function add(\n    f1     :: MatsubaraFunction{GD, SD, DD, Q}, \n    f2     :: MatsubaraFunction{GD, SD, DD, Q}\n    ;\n    checks :: Bool = true\n    )      :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}\n\nAddition of two MatsubaraFunction, returns new MatsubaraFunction. Safety measures can be disabled with checks = false if needed for performance (discouraged).\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"add!","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.add!","page":"MatsubaraFunction","title":"MatsubaraFunctions.add!","text":"function add!(\n    f1     :: MatsubaraFunction{GD, SD, DD, Q}, \n    f2     :: MatsubaraFunction{GD, SD, DD, Q}\n    ;\n    checks :: Bool = true\n    )      :: Nothing where {GD, SD, DD, Q <: Number}\n\nInplace addition of two MatsubaraFunction (f1 += f2). Safety measures can be disabled with checks = false if needed for performance (discouraged).\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"subtract","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.subtract","page":"MatsubaraFunction","title":"MatsubaraFunctions.subtract","text":"function subtract(\n    f1     :: MatsubaraFunction{GD, SD, DD, Q}, \n    f2     :: MatsubaraFunction{GD, SD, DD, Q}\n    ;\n    checks :: Bool = true\n    )      :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}\n\nSubtraction of two MatsubaraFunction, returns new MatsubaraFunction. Safety measures can be disabled with checks = false if needed for performance (discouraged).\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"subtract!","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.subtract!","page":"MatsubaraFunction","title":"MatsubaraFunctions.subtract!","text":"function subtract!(\n    f1     :: MatsubaraFunction{GD, SD, DD, Q}, \n    f2     :: MatsubaraFunction{GD, SD, DD, Q}\n    ;\n    checks :: Bool = true\n    )      :: Nothing where {GD, SD, DD, Q <: Number}\n\nInplace subtraction of two MatsubaraFunction (f1 -= f2). Safety measures can be disabled with checks = false if needed for performance (discouraged).\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"mult","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.mult","page":"MatsubaraFunction","title":"MatsubaraFunctions.mult","text":"function mult(\n    f      :: MatsubaraFunction{GD, SD, DD, Q},\n    val    :: Qp\n    ;\n    checks :: Bool = true\n    )      :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number, Qp <: Number}\n\nMultiplication of MatsubaraFunction with scalar, returns new MatsubaraFunction. Safety measures can be disabled with checks = false if needed for performance (discouraged).\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"mult!","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.mult!","page":"MatsubaraFunction","title":"MatsubaraFunctions.mult!","text":"function mult!(\n    f   :: MatsubaraFunction{GD, SD, DD, Q},\n    val :: Qp\n    )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}\n\nInplace multiplication of MatsubaraFunction with scalar (f *= val)\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"set!","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.set!","page":"MatsubaraFunction","title":"MatsubaraFunctions.set!","text":"function set!(\n    f   :: MatsubaraFunction{GD, SD, DD, Q},\n    val :: Qp,\n    )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}\n\nInitialize MatsubaraFunction with val\n\n\n\n\n\nfunction set!(\n    f1     :: MatsubaraFunction{GD, SD, DD, Q},\n    f2     :: MatsubaraFunction{GD, SD, DD, Q},\n    ; \n    checks :: Bool = true\n    )      :: Nothing where {GD, SD, DD, Q <: Number}\n\nInitialize MatsubaraFunction with another MatsubaraFunction (f1 = f2). Safety measures can be disabled with checks = false if needed for performance (discouraged).\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"LinearIndex","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.LinearIndex","page":"MatsubaraFunction","title":"MatsubaraFunctions.LinearIndex","text":"function LinearIndex(\n    f :: MatsubaraFunction{GD, SD, DD, Q},\n    w :: NTuple{GD, MatsubaraFrequency},\n    x :: Vararg{Int64, SD} \n    ) :: Int64 where {GD, SD, DD, Q <: Number}\n\nReturns linear index for access to f.data\n\n\n\n\n\nfunction LinearIndex(\n    f    :: MatsubaraFunction{GD, SD, DD, Q},\n    cidx :: CartesianIndex{DD}\n    )    :: Int64 where {GD, SD, DD, Q <: Number}\n\nReturns linear index for access to f.data\n\n\n\n\n\nfunction LinearIndex(\n    f :: MatsubaraFunction{GD, SD, DD, Q},\n    x :: Vararg{Int64, DD}\n    ) :: Int64 where {GD, SD, DD, Q <: Number}\n\nReturns linear index for access to f.data    \n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"to_Matsubara","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.to_Matsubara","page":"MatsubaraFunction","title":"MatsubaraFunctions.to_Matsubara","text":"function to_Matsubara(\n    f    :: MatsubaraFunction{GD, SD, DD, Q},\n    cidx :: CartesianIndex{DD}\n    )    :: Tuple{NTuple{GD, MatsubaraFrequency}, NTuple{SD, Int64}} where {GD, SD, DD, Q <: Number}\n\nReturns coordinates in grids and index of tensor structure\n\n\n\n\n\nfunction to_Matsubara(\n    f   :: MatsubaraFunction{GD, SD, DD, Q},\n    idx :: Int64 \n    )   :: Tuple{NTuple{GD, MatsubaraFrequency}, NTuple{SD, Int64}} where {GD, SD, DD, Q <: Number}\n\nReturns coordinates in grids and index of tensor structure\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"upper_tail_moments","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.upper_tail_moments","page":"MatsubaraFunction","title":"MatsubaraFunctions.upper_tail_moments","text":"function upper_tail_moments(\n    f  :: MatsubaraFunction{1, SD, DD, Q},\n    α0 :: Q,\n    x  :: Vararg{Int64, SD} \n    )  :: Tuple{Q, Q} where {SD, DD, Q <: Number}\n\nReturns high frequency moments for quadratic model using upper grid bound. Here, α0 is the asymptotic limit for large positive frequencies. Note, that the distance between interpolation nodes  is kept constant below T = 1.\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"lower_tail_moments","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.lower_tail_moments","page":"MatsubaraFunction","title":"MatsubaraFunctions.lower_tail_moments","text":"function lower_tail_moments(\n    f  :: MatsubaraFunction{1, SD, DD, Q},\n    α0 :: Q,\n    x  :: Vararg{Int64, SD} \n    )  :: Tuple{Q, Q} where {SD, DD, Q <: Number}\n\nReturns high frequency moments for quadratic model using lower grid bound. Here, α0 is the asymptotic limit for large negative frequencies. Note, that the distance between interpolation nodes  is kept constant below T = 1.\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"sum_me","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.sum_me","page":"MatsubaraFunction","title":"MatsubaraFunctions.sum_me","text":"function sum_me(\n    f  :: MatsubaraFunction{1, SD, DD, Q},\n    α0 :: Q,\n    x  :: Vararg{Int64, SD}\n    )  :: Q where {SD, DD, Q <: Complex}\n\nComputes the Matsubara sum (with regulator exp(-iw0+)) for a complex valued MatsubaraFunction on 1D grid. Here, α0  is the asymptotic limit for large frequencies. This is only viable if f has a Laurent series representation with respect  to an annulus about the imaginary axis.\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"sgn","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.sgn","page":"MatsubaraFunction","title":"MatsubaraFunctions.sgn","text":"sgn(op :: MatsubaraOperation) :: Bool\n\nReturn op.sgn\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"con","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.con","page":"MatsubaraFunction","title":"MatsubaraFunctions.con","text":"con(op :: MatsubaraOperation) :: Bool\n\nReturn op.con\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"save_matsubara_function!","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.save_matsubara_function!","page":"MatsubaraFunction","title":"MatsubaraFunctions.save_matsubara_function!","text":"function save_matsubara_function!(\n    h :: HDF5.File,\n    l :: String,\n    f :: MatsubaraFunction{Dg, Ds, Dt, Q}\n    ) :: Nothing where {Dg, Ds, Dt, Q <: Number}\n\nSave MatsubaraFunction f with label l to file h   \n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"load_matsubara_function","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.load_matsubara_function","page":"MatsubaraFunction","title":"MatsubaraFunctions.load_matsubara_function","text":"function load_matsubara_function(\n    h :: HDF5.File,\n    l :: String\n    ) :: MatsubaraFunction\n\nLoad MatsubaraFunction with label l from file h\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"save_matsubara_symmetry_group!","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.save_matsubara_symmetry_group!","page":"MatsubaraFunction","title":"MatsubaraFunctions.save_matsubara_symmetry_group!","text":"function savematsubarasymmetry_group!(     h  :: HDF5.File,     l  :: String,     SG :: MatsubaraSymmetryGroup     )  :: Nothing\n\nSave MatsubaraSymmetryGroup SG with label l to file h\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"load_matsubara_symmetry_group","category":"page"},{"location":"matsubara_func/#MatsubaraFunctions.load_matsubara_symmetry_group","page":"MatsubaraFunction","title":"MatsubaraFunctions.load_matsubara_symmetry_group","text":"function loadmatsubarasymmetry_group(     h :: HDF5.File,     l :: String     ) :: MatsubaraSymmetryGroup\n\nLoad MatsubaraSymmetryGroup with label l from file h\n\n\n\n\n\n","category":"function"},{"location":"matsubara_func/","page":"MatsubaraFunction","title":"MatsubaraFunction","text":"","category":"page"},{"location":"#MatsubaraFunctions.jl","page":"Home","title":"MatsubaraFunctions.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package aims at providing a convenient interface to rapidly prototype algorithms for multivariable Green's functions of the form G_i_1  i_n(iomega_1  iomega_m), where i_k denote lattice or orbital indices and omega_l are fermionic/bosonic Matsubara frequencies.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package is not yet registered, but available from","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add https://github.com/dominikkiese/MatsubaraFunctions.jl","category":"page"},{"location":"io/#IO","page":"IO","title":"IO","text":"","category":"section"}]
}
