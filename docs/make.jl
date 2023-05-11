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
            "IO"                 => "io.md"])

deploydocs(; repo = "github.com/dominikkiese/MatsubaraFunctions.jl",)