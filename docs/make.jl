push!(LOAD_PATH, "../src/")
import Pkg; Pkg.add("Documenter")
Pkg.add("HDF5")

using MatsubaraFunctions
using Documenter
using HDF5

makedocs(sitename = "MatsubaraFunctions.jl", 
          modules = [MatsubaraFunctions],
          pages   = [
            "Home"         => "index.md", 
            "MeshFunction" => "mesh_func.md",

            "Mesh" => ["mesh.md",
              "Specific meshes" => [
                "matsubara_mesh.md",
                "brillouin_zone_mesh.md",
                "index_mesh.md",
                "custom_mesh.md"
              ]],

            "Advanced usage" => [
              "Symmetry reduction" => "mesh_func_symmetries.md",
              "Parallelization"    => "mesh_func_parallelization.md"
            ],

            "Miscellaneous" => "misc.md"])


deploydocs(; repo = "github.com/dominikkiese/MatsubaraFunctions.jl",)