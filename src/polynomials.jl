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
    legendre_recurrence(x, n)

Compute legendre polynomial at point `x` of order `n_plus_1` using 
recurrence relation and dynamic programming.

# Examples
```julia-repl
julia> using PDEMethods: legendre_recurrence
julia> using Plots
julia> x = 2
julia> @assert legendre_recurrence(x, 0) == 1
julia> @assert legendre_recurrence(x, 1) == x
julia> @assert legendre_recurrence(x, 2) == 0.5*(3*x^2 - 1)
julia> @assert legendre_recurrence(x, 3) == 0.5*(5*x^3 - 3*x)
julia> x_range = -1:0.01:1
julia> myp = plot(
    x_range, legendre_recurrence.(x_range, 0), label = "P_0", ylims = (-1, 1));
julia> for order = 1:5
            plot!(
                myp, x_range, legendre_recurrence.(x_range, order), 
                label="P_" * string(order))
       end
julia> myp
```

# References
[1] Sullivan2015 pg. 140 "8.2: Recurrence Relations"
"""
function legendre_recurrence(x, N)
    L = zeros(N+1)
    L[1] = 1           # L[1] = L_0
    if N >= 1 
        L[2] = x       # L[2] = L_1
        for i = 3:N+1  # L[3] = L_2
            ix = i - 1 # Must use previously stored evaluations of L     
            order_n = ix - 1  # Order L_n is the ix minus 1
            L[i] = ((2*(order_n)+1)*x*L[ix] - order_n*L[ix-1])/(order_n+1)
        end
    end 
    return L[end]
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
        interp_polynomial += y[j]*lagrange_basis(p, j, x)
    end
    return interp_polynomial
end 

"""
    lagrange_basis(p, j, x)

Return lagrange basis polynomial at point `p` about a nodal index 
`j` for a set of discrete x-values (`x`).
"""
function lagrange_basis(p, j, x)
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
