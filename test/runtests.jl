#using Revise
using Test
using Aqua
using HDF5
using MPI
using StaticArrays
using OffsetArrays
using MatsubaraFunctions 

MPI.Init()
MatsubaraFunctions.DEBUG() = true # enable all checks for testing
Aqua.test_all(MatsubaraFunctions)

test_dir = joinpath(dirname(@__FILE__), "tests")
for test in readdir(test_dir)
    include(joinpath(test_dir, test))
end