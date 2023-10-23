"""
    legendre(x, n)
  
Compute Legendre polynomial of order `n` at point `x` using explicit
representation.

# References
[1] : https://en.wikipedia.org/wiki/Legendre_polynomials
"""
function legendre(x, n)
    @assert n >= 0
    polysum = 0
    for k = 0:n
        intermed = binomial(n, k)*binomial(n+k, k)*((x - 1)/2)^k
        polysum += intermed
    end 
    return polysum 
end

"""
    hermite(x, n)

Compute physicist's Hermite polynomial of order `n` at point `x` using 
explicit expression.

# References
[1] : https://en.wikipedia.org/wiki/Hermite_polynomials
"""
function hermite(x, n)
    @assert n >= 0
    polysum = 0
    upper_bound = Int(floor(n/2))
    for m = 0:upper_bound
        polysum += ((2*x)^(n-2*m)*(-1)^m)/(factorial(m)*factorial(n - 2*m))
    end 
    return factorial(n)*polysum
end

"""
    lagrange(p, x, y)

Return lagrange interpolating polynomial.

# Examples
```julia-repl
julia> using Plots
julia> using PDEMethods: lagrange
julia> x = [-9, -4, -1, 7]
julia> y = [5, 2, -2, 9]
julia> pl = scatter(x, y, color = :red, legend = false)
julia> continuous_x = LinRange(x[1], x[end], 100)
julia> lagrange_ys = lagrange.(continuous_x, Ref(x), Ref(y))
julia> plot!(pl, continuous_x, lagrange_ys, legend = false)
```

# References
[1] : https://en.wikipedia.org/wiki/Lagrange_polynomial
"""
function lagrange(p, x, y) 
    k_plus_one = length(x)
    @assert k_plus_one == length(y) "`x` and `y` are coordinate pairs"
    interp_polynomial = 0
    for j = 1:k_plus_one
        interp_polynomial += y[j]*lagrange_basis(p, j, x, y)
    end
    return interp_polynomial
end 

"""
    lagrange_basis(p, j, x, y)

Return lagrange basis polynomial at point `p` about a nodal index 
`j` for a set of discrete x-values (`x`) and corresponding y-values (`y`).
"""
function lagrange_basis(p, j, x, y)
    #   1   2  ...   k+1
    # {x0, x1, ... , xk} |> length --> k + 1
    k_plus_one = length(x)

    # basis polynomial product
    basis_polynomial = 1
    for m = 1:k_plus_one
        m != j && (basis_polynomial *= (p - x[m])/(x[j] - x[m]))
    end

    return basis_polynomial
end
