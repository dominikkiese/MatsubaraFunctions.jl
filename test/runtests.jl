using Test
using MPI
using HDF5
using MatsubaraFunctions 
using Aqua

MPI.Init()
MatsubaraFunctions.DEBUG() = true # enable all checks for testing
Aqua.test_all(MatsubaraFunctions)

for file in readdir("tests")
    include("tests/" * file)
end