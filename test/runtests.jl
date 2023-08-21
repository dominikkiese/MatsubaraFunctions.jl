using Test
using MPI; MPI.Init()
using HDF5
using MatsubaraFunctions 

for file in readdir("tests")
    include("tests/" * file)
end