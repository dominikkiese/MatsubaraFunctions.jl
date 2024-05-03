push!(LOAD_PATH, "../src/")
import Pkg; Pkg.add("Documenter")

using MatsubaraFunctions
using Documenter

makedocs(sitename = "MatsubaraFunctions.jl", 
          modules = [MatsubaraFunctions],
          pages   = [
            "Home"               => "index.md", 
            "MeshFunction"  => "mesh_func.md",
            "Mesh"      => "mesh.md",
            "MeshPoints" => [
              "brillouin_zone_mesh.md",
              "index_mesh.md",
              "matsubara_mesh.md",
            ],
            "Advanced usage"     => [
              "Symmetry reduction"   => "mesh_func_symmetries.md",
              "Parallelization"      => "mesh_func_parallelization.md",
              #"Matsubara sums"       => "matsubara_func_sums.md"
            ],
            "Miscellaneous"      => "misc.md"])

deploydocs(; repo = "github.com/dominikkiese/MatsubaraFunctions.jl",)