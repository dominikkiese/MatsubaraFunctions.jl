push!(LOAD_PATH, "../src/")
import Pkg; Pkg.add("Documenter")

using MatsubaraFunctions
using Documenter

makedocs(sitename = "MatsubaraFunctions.jl", 
          modules = [MatsubaraFunctions],
          pages   = ["Home" => "index.md"])

deploydocs(; repo = "github.com/dominikkiese/MatsubaraFunctions.jl",)