# Iterative methods for solving linear systems begin with an initial estimate 
# for the solution and successively improve it until the solution is as accurate 
# as desired:
# * jacobi
# * gauss-seidel
# * conjugate gradient method (+ preconditioning?)

"""
    jacobi_method(A::Matrix, b::Vector, x0::Vector, niters::Int)

```math
x_{i}^{(k+1)} = \\frac{b_i - \\sum_{j \\neq i} a_{ij} x_j^{(k)}}{a_{ii}}, i \\in [1...n] 
```

Return iterative solution `x` to linear system `Ax = b` for coefficient matrix 
`A`, rhs vector `b`, and initial guess `x0` using Jacobi method for `niters`.

# Examples
```julia-repl
julia> using PDEMethods: jacobi_method, laplace_eq_init_arrays
julia> mesh, A, b = laplace_eq_init_arrays(grid_dim_n=4, lbc=0, rbc=0, tbc=1, bbc=0)
julia> direct_u = A \\b
julia> jacobi_u = jacobi_method(A, b, zeros(length(b)), 100)
julia> @assert all(direct_u .≈ jacobi_u) 
```

# References
[1] : Heath ch. 11.5.2
"""
function jacobi_method(A::Matrix, b::Vector, x0::Vector, niters::Int)
    xkplus1 = zeros(length(b))
    xk = x0
    nrows, ncols = size(A)
    for k in 1:niters
        for i in 1:nrows
            prod_sum = 0
            for j in 1:ncols
                i != j && (prod_sum += A[i, j]*xk[j])
            end
            xkplus1[i] = (b[i] - prod_sum)/A[i, i]
        end
        xk = xkplus1
    end
    return xkplus1
end 
