# Advanced Usage: Matsubara Sums 

!!! warning "sum_me deprecated"
    The function sum_me has been deprecated and will not be supported in future releases.

For `MatsubaraFunction` objects $G_{i_1 ... i_n}(i\omega)$ defined on 1D grids, we export the function `sum_me`, which computes the series $\Sigma_m G_{i_1 ... i_n}(i\omega_{m}) e^{i\omega_m 0^+}$ for $m \in \mathbb{Z}$ using tail fits of $G$ together with analytic formulas for summations of the form $\Sigma_m \frac{1}{(i\omega_m)^\alpha}e^{i\omega_m 0^+}$ with $\alpha \in \mathbb{N}$. This, however, requires $G$ to be representable by a Laurent series in an elongated annulus about the imaginary axis. Also, $G$ must decay to zero. This feature is experimental and the API may change in future versions.

```julia
ξ = 0.5
T = 1.0
N = 128
g = MatsubaraGrid(T, N, Fermion)
f = MatsubaraFunction(g, 1)

for v in g
    f[v] = 1.0 / (im * value(v) - ξ)
end 

# evaluate the series and compare to analytic result
ρ(x, T) = 1.0 / (exp(x / T) + 1.0)
println(abs(sum_me(f) - (ρ(+ξ, T) - 1.0)))
```