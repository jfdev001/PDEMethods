# Direct methods for sparse linear systems 
# should annihilate sparse entries (minimum degree algorithms)
# * LU (Gaussian Elimination)
# * Cholesky Factorization

using LinearAlgebra: UnitLowerTriangular, UpperTriangular, LowerTriangular,
    issymmetric, isposdef

"""
    lu_factorization!(A)

Return LU factorization of matrix `A` into lower triangular matrix `L` and upper 
triangular matrix `U`. If `inplace = true`, then the subdiagonal entries
of `A` are the lower triangular matrix `L` while the diagonal and superdiagonal
entries of `A` are the upper triangular matrix `U`.

# References 
[1] : Heath algorithm 2.3
"""
function lu_factorization!(A)
    m, n = size(A)
    for k = 1:n-1
        A[k, k] == 0 && break
        for i = k+1:n
            multiplier =  A[i, k]/A[k, k]
            A[i, k] = multiplier 
        end 

        for j = k+1:n
            for i = k+1:n
                submatrix_transform = A[i, j] - A[i, k]*A[k, j] 
                A[i, j] = submatrix_transform 
            end
        end
    end
    return nothing
end 

"""
    lu_factorization(A)

Return `L` and `U` factorized matrices of `A` explicitly via OOP operation.
"""
function lu_factorization(A) 
    A_ = copy(A)
    lu_factorization!(A_)
    # Explicitly cast to matrices to show zeros
    L = Matrix(UnitLowerTriangular(A_))
    U = Matrix(UpperTriangular(A_))
    return L, U
end

"""
    cholesky_factorization!(A)

Return `L` such that `A = L*transpose(L)` where `L` is lower triangular
but does not, in general, have a unit diagonal.

# References
[1] : Heath Algorithm 2.7
"""
function cholesky_factorization!(A)
    @assert issymmetric(A) && isposdef(A) "`A` must be SPD"
    m, n = size(A)
    for k = 1:n
        A[k, k] = sqrt(A[k, k])
        for i = k+1:n
            A[i, k] = A[i, k]/A[k, k]
        end 

        for j = k+1:n
            for i = j:n
                A[i, j] = A[i, j] - A[i, k]*A[j, k]
            end
        end
    end
    
    # (inefficient) set superdiagonal entries to 0
    for j = 2:n
        for i = 1:m
            if i < j
                A[i, j] = 0.0
            end
        end
    end 
    return nothing
end

"""
    cholesky_factorization(A)

Return `L` from OOP cholesky factorization.

# Examples
```julia-repl
julia> using PDEMethods: cholesky_factorization, laplace_eq_init_arrays
julia> using PDEMethods: forward_substitution, back_substitution
julia> mesh, A, b = laplace_eq_init_arrays()
julia> L = cholesky_factorization(A) 
julia> Lᵀ = transpose(L)              # LLᵀx = b
julia> y = forward_substitution(L, b) # L y = b === L \\ b
julia> x = back_substitution(Lᵀ, y)   # Lᵀx = y === Lᵀ \\ y 
julia> @assert all(x .≈ A \\ b )
```
"""
function cholesky_factorization(A) 
    L = copy(A)
    cholesky_factorization!(L)
    return L
end

"""   
    forward_substitution(L::Matrix, b_::Vector)

Solve a system of equations `Ly = b` with lower triangular system `L` and 
rhs `b` by forward substitution.
 
# References
[1] : Heath algorithm 2.1.
"""
function forward_substitution(L::AbstractMatrix, b_::Vector)
    # L is square, so m, n doesn't really matter here
    m, n = size(L)
    y = zeros(n)
    b = copy(b_)
    for j = 1:n # loops over columns
        if L[j, j] == 0
            @show break
        end
        y[j] = b[j]/L[j, j]
        for i = (j + 1):n
            b[i] = b[i] - L[i, j]*y[j]
        end
    end 
    return y
end


"""
    back_substitution(U::Matrix, y_::Vector)

Solve a system `Ux = y` by backward substitution of the upper triangular
system `U` and rhs vector `y_`.

# References
[1] : Heath algorithm 2.2, replaced `b` with `y_`.
"""
function back_substitution(U::AbstractMatrix, y_::Vector)
    m, n = size(U)
    y = copy(y_)
    x = zeros(n)
    for j = n:-1:1
        U[j,j] == 0 && break
        x[j] = y[j]/U[j, j]
        for i = 1:j-1
            y[i] = y[i] - U[i, j]*x[j]
        end
    end
    return x
end


function back_substitution_matheq(U::AbstractMatrix, b::Vector) 
    n, m = size(U) 
    x = zeros(n)
    x[n] = b[n]/U[n, n]
    for i = n-1:-1:1
        inner_sum = 0
        for j = i+1:n
            inner_sum += U[i, j]*x[j]
        end
        U[i, i] == 0 && break
        x[i] = (b[i] - inner_sum)/U[i,i]
    end
    return x
end 


"""
    forward_substitution_matheq(L::Matrix, b::Vector)

Direct conversion of equation in Heath pg. 64 describing forward substitution.
"""
function forward_substitution_matheq(L::AbstractMatrix, b::Vector)
    n, m = size(L) # n == nrows based on eq, tho n == m because square anyway
    x = zeros(n) 
    x[1] = b[1]/L[1, 1]
    for i = 2:n
        inner_sum = 0
        for j = 1:i-1
            inner_sum += L[i, j]*x[j]
        end
        if L[i, i] == 0 # prevent divide 0
            @show break
        end
        (x[i] = (b[i] - inner_sum)/L[i, i])
    end
    return x
end

