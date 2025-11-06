# Go through the publication of MatsubaraFunctions.jl a write a version of the code examples compatible with the unstable branch.

using MatsubaraFunctions

# Section 3.1
T = 1.0
k = 10
v = MatsubaraFrequency(T, k, Fermion)
w = MatsubaraFrequency(T, k, Boson)