# Iterative methods for solving linear systems begin with an initial estimate 
# for the solution and successively improve it until the solution is as accurate 
# as desired:
# * jacobi
# * jacobi based on block diagonal formula only (see Dolean2015)
# * block jacobi 
# * gauss-seidel
# * conjugate gradient method (+ preconditioning?)

using LinearAlgebra: LowerTriangular, UpperTriangular, diag, diagm, diagind,
    Diagonal

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
                i != j && (prod_sum += A[i, j] * xk[j])
            end
            xkplus1[i] = (b[i] - prod_sum) / A[i, i]
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
[2] : Ch.1.2 Dolean2015
"""
function jacobi_method_dolean(A::Matrix, b::Vector, x0::Vector, niters::Int)
    xk = x0
    #D = block_diagonal_matrix(A)
    D_vec = diag(A)
    D = diagm(D_vec)
    D_inv = inv(D)
    for k in 1:niters
        rk = b - A * xk
        xk = xk + D_inv * rk
    end
    return xk
end

function jacobi_method_matrix_form_solver(
    A::Matrix, b::Vector, x0::Vector, niters::Int)
    x = x0
    D_vec = diag(A)
    D = diagm(D_vec)
    D_inv = inv(D)
    L = copy(LowerTriangular(A)) # constructs a view, copy to avoid change A
    L[diagind(L)] .= 0
    U = copy(UpperTriangular(A)) # constructs a view, copy to avoid change A
    U[diagind(U)] .= 0
    for k in 1:niters
        x = D_inv * (b - (L + U) * x)
    end
    return x
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
function block_diagonal_matrix(A::Matrix, block_size::Int=2)
    m, n = size(A)
    @assert m == n "square matrix"
    M = zeros(m, n)
    diagonal_col_ix = 1
    for i in 1:m
        diagonal_col_ix += (i != 1 && isodd(i)) * block_size
        for j in diagonal_col_ix:diagonal_col_ix+(block_size-1)
            M[i, j] = A[i, j]
        end
    end
    return M
end

"""
    gauss_seidel_matheq(A::Matrix, b::Vector, x0::Vector, niters::Int; ω = 1.0)

Return solution to linear system `Ax = b` via Gauss-Seidel (ω = 1., or SOR
for ω != 1 and ω ∈ (0, 2)) method.

# References
[1] : Ch. 11.5.3 Heath.
"""
function gauss_seidel_matheq(
    A::Matrix, b::Vector, x0::Vector, niters::Int; ω=1.0)
    @assert 0 < ω < 2
    x = x0
    m, n = size(A)
    for k in 1:niters
        for i in 1:n
            xkplus1_sum = 0
            for j in 1:i-1
                xkplus1_sum += A[i, j] * x[j]
            end

            xk_sum = 0
            for j in i+1:n
                xk_sum += A[i, j] * x[j]
            end
           
            xi_kplus1_GS = (b[i] - xkplus1_sum - xk_sum) / A[i, i] # Gauss Seidel method 
            xi_k = x[i]
            x[i] = xi_k + ω*(xi_kplus1_GS - xi_k) # update x_i^{(k+1)}
        end
    end
    return x
end


"""
    gauss_seidel_matrix_form_solver(
        A::Matrix, b::Vector, x0::Vector, niters::Int)

Return solution to linear system `Ax = b` via matrix Gauss-Seidel method.

# References
[1] : Ch. 11.5.3 Heath.
"""
function gauss_seidel_matrix_form_solver(
    A::Matrix, b::Vector, x0::Vector, niters::Int)
    # Get diagonal entries of A, store in matrix
    D_vec = diag(A)
    D = diagm(D_vec)

    # Get the strict lower and upper triangular matrices
    L = LowerTriangular(A)
    L[diagind(L)] .= 0
    U = UpperTriangular(A)
    U[diagind(U)] .= 0

    # Compute the inverse of D + L ahead of time
    D_plus_L_inv = inv(D + L)

    # initialize the guess
    x = x0

    # gauss seidel iterations
    for k in 1:niters
        x = D_plus_L_inv * (b - U * x)
    end

    return x
end

function ssor_preconditioner(A_::Matrix, ω::Float64=1.0)
    @assert issymmetric(A_)
    @assert 0 < ω < 2
    A = copy(A_) # expensive, whatevs
    D = Diagonal(A)
    L = LowerTriangular(A)
    Lᵀ = transpose(L)
    M = (1 / (2 - ω)) * ((1 / ω) * D + L) * inv(((1 / ω) * D)) * transpose((1 / ω) * D + L)
    return M
end

"""
    cg_poisson(nx, ny, un, r, max_iter, init_rms = nothing)

Conjugate gradient method for 2D Poisson equation. 

# Arguments
- `nx, ny`: number of grid points in x and y directions
- `nx, ny`: grid spacing in x and y directions
- `un`: numerical solution matrix for the Poisson equation
- `r`: residual matrix
- `max_iter`: maximum number if iterations
- `init_rms`: L-2 norm of the residual at the start of iterations

# References
[1] : Listing 22 in CFD Julia: CFD Julia: A Learning Module Structuring an 
Introductory Course on Computational Fluid Dynamics
"""
function cg_poisson(nx, ny, un, r, max_iter, initial_rms)
    p = r # assign initial residual to the conjugate vector
    for k = 1:max_iter

        # conjugate gradient algorithm
        for j = 2:ny
            for i = 2:nx
                # centered difference formulas for second order deriv in x and y dir
                q[i, j] = (p[i+1, j] - 2.0 * p[i, j] + p[i-1, j]) / (dx^2) +
                        (p[i, j+1] - 2.0 * p[i, j] + p[i, j-1]) / (dy^2)
            end
        end

        aa, bb = 0.0, 0.0
        for j = 2:ny
            for i = 2:nx
                aa = aa + r[i, j] * r[i, j] # <r,r>
                bb = bb + q[i, j] * p[i, j] # <q,p>
            end
        end

        cc = aa / (bb + tiny) #alpha
        for j = 2:ny
            for i = 2:nx
                un[i, j] = un[i, j] + cc * p[i, j] # update the numerical solution
            end
        end

        bb, aa = aa, 0.0
        for j = 2:ny
            for i = 2:nx
                r[i, j] = r[i, j] - cc * q[i, j] # update the residual
                aa = aa + r[i, j] * r[i, j]
            end
        end

        cc = aa / (bb + tiny) # beta
        for j = 1:ny
            for i = 1:nx
                p[i, j] = r[i, j] + cc * p[i, j] # update the conjugate vector
            end
        end

        # compute the residual at k-th step using the new solution
        for j = 2:ny
            for i = 2:nx
                d2udx2 = (un[i+1, j] - 2 * un[i, j] + un[i-1, j]) / (dx^2)
                d2udy2 = (un[i, j+1] - 2 * un[i, j] + un[i, j-1]) / (dy^2)
                r[i, j] = f[i, j] - d2udx2 - d2udy2
            end
        end

        # calculate the L-2 norm of the residual vector
        rms = 0.0
        for j = 2:ny
            for i = 2:nx
                rms = rms + r[i, j]^2
            end
        end

        rms = sqrt(rms / ((nx - 1) * (ny - 1)))
        # if the convergence critera (rms<tolerance) is satidfied, stop the iteration
        if (rms / initial_rms) <= 1e-12
            break
        end
    end
end 
