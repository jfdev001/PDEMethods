# Iterative methods for solving linear systems begin with an initial estimate 
# for the solution and successively improve it until the solution is as accurate 
# as desired:
# * jacobi
# * jacobi based on block diagonal formula only (see Dolean2015)
# * block jacobi 
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

"""
    jacobi_method_dolean(A::Matrix, b::Vector, x0::Vector, niters::Int)

Return solution to linear system of equation.

# Examples
```julia-repl
julia> using PDEMethods: jacobi_method_dolean, laplace_eq_init_arrays
julia> mesh, A, b = laplace_eq_init_arrays(grid_dim_n=4, lbc=0, rbc=0, tbc=1, bbc=0)
julia> direct_u = A \\ b
julia> jacobi_u = jacobi_method_dolean(A, b, zeros(length(b)), 100)
julia> @assert all(direct_u .≈ jacobi_u) 
```

# References
[1] : [The Jacobi Algorithm](https://www.youtube.com/watch?v=hcoikfp64bM&list=TLPQMTQxMTIwMjORT8jzlo-9xw&index=1)
"""
function jacobi_method_dolean(A::Matrix, b::Vector, x0::Vector, niters::Int)
    xk = x0
    D = block_diagonal_matrix(A)
    D_inv = inv(D)
    for k in 1:niters
        rk = b - A*xk
        xk = xk + D_inv*rk
    end  
    return xk
end

function block_jacobi_method(
    A::Matrix, p_procs::Int, b::Vector, x0::Vector, niters::Int)

    throw("notimplementederror")
    m, n = size(A)
    @assert m >= p_procs
    remaining_rows = mod(m, p_procs)
    
    restriction_matrices = nothing
end 

"""
    block_diagonal_matrix(A::Matrix, block_size::Int = 2)

Return the block diagonal matrix of the matrix `A`.

# Examples
```julia-repl
julia> using PDEMethods: laplace_eq_init_arrays, block_diagonal_matrix
julia> mesh, A, b = laplace_eq_init_arrays(;grid_dim_n=4, lbc=0, rbc=0, tbc=1, bbc=0)
julia> A
4×4 Matrix{Float64}:
  4.0  -1.0  -1.0   0.0
 -1.0   4.0   0.0  -1.0
 -1.0   0.0   4.0  -1.0
  0.0  -1.0  -1.0   4.0
julia> block_diagonal_matrix(A) 
4×4 Matrix{Float64}:
  4.0  -1.0   0.0   0.0
 -1.0   4.0   0.0   0.0
  0.0   0.0   4.0  -1.0
  0.0   0.0  -1.0   4.0
```

# References
[1] : https://phtournier.pages.math.cnrs.fr/5mm29/blockjacobi/
"""
function block_diagonal_matrix(A::Matrix, block_size::Int = 2)
    m, n = size(A)
    @assert m == n "square matrix"
    M = zeros(m, n)
    diagonal_col_ix = 1
    for i in 1:m
        diagonal_col_ix += (i != 1 && isodd(i))*block_size
        for j in diagonal_col_ix:diagonal_col_ix+(block_size-1)
           M[i, j] = A[i, j]
        end
    end
    return M 
end 

"""
    gauss_seidel_matheq(A::Matrix, x0::Vector, b::Vector, niters::Int)

Return solution to linear system `Ax = b` via Gauss-Seidel method.

# References
[1] : Ch. 11.5.3 Heath.
"""
function gauss_seidel_matheq(A::Matrix, x0::Vector, b::Vector, niters::Int)
    x = x0
    m, n = size(A)
    for k in 1:niters
        for i in 1:n
            xkplus1_sum = 0
            for j in 1:i-1
                xkplus1_sum += A[i, j]*x[j]
            end
        
            xk_sum = 0
            for j in i+1:n
                xk_sum += A[i, j]*x[j]
            end 

            x[i] = (b[i] - xkplus1_sum - xk_sum)/A[i,i]
        end
    end
    return x
end 
