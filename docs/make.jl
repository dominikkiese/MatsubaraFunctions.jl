push!(LOAD_PATH, "../src/")
using MatsubaraFunctions
import Pkg; Pkg.add("Documenter")
using Documenter

makedocs(sitename = "MatsubaraFunctions.jl", 
          modules = [MatsubaraFunctions],
          pages   = ["Home" => "index.md"])

deploydocs(; repo = "github.com/dominikkiese/MatsubaraFunctions.jl",)