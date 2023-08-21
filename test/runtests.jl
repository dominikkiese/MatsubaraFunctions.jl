using Test
using MPI; MPI.Init()
using HDF5
using MatsubaraFunctions 
using Aqua

Aqua.test_all(MatsubaraFunctions)

for file in readdir("tests")
    include("tests/" * file)
end