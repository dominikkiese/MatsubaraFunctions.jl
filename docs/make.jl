push!(LOAD_PATH, "../src/")
import Pkg; Pkg.add("Documenter")

using MatsubaraFunctions
using Documenter

makedocs(sitename = "MatsubaraFunctions.jl", 
          modules = [MatsubaraFunctions],
          pages   = [
            "Home"               => "index.md", 
            "MatsubaraFrequency" => "matsubara_freq.md",
            "MatsubaraGrid"      => "matsubara_grid.md",
            "MatsubaraFunction"  => "matsubara_func.md",
            "Advanced usage"     => [
              "Symmetry reduction"   => "matsubara_func_symmetries.md",
              "Parallelization"      => "matsubara_func_parallelization.md",
              "Matsubara sums"       => "matsubara_func_sums.md"
            ],
            "Miscellaneous"      => "misc.md",
            "MatsubaraIndex"     => "matsubara_index.md"])

deploydocs(; repo = "github.com/dominikkiese/MatsubaraFunctions.jl",)