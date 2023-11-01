# Direct methods for sparse linear systems 
# should annihilate sparse entries (minimum degree algorithms)
# * LU (Gaussian Elimination)
# * Cholesky Factorization

using LinearAlgebra: UnitLowerTriangular, UpperTriangular, LowerTriangular

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
"""
function cholesky_factorization(A) 
    L = copy(A)
    cholesky_factorization!(L)
    return L
end

function forward_substitution

end
