using Test
using Aqua
using HDF5
using MPI
using StaticArrays
using MatsubaraFunctions 

MPI.Init()
MatsubaraFunctions.DEBUG() = true # enable all checks for testing
Aqua.test_all(MatsubaraFunctions)

for test in readdir("tests")
    include(joinpath("tests", test))
end