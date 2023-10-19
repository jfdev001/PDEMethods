# Direct methods for sparse linear systems 
# should annihilate sparse entries (minimum degree algorithms)
# * LU (Gaussian Elimination)
# * Cholesky Factorization

using LinearAlgebra: UnitLowerTriangular, UpperTriangular

"""
    gaussian_elimination!(A)

Return LU factorization of matrix `A` into lower triangular matrix `L` and upper 
triangular matrix `U`. If `inplace = true`, then the subdiagonal entries
of `A` are the lower triangular matrix `L` while the diagonal and superdiagonal
entries of `A` are the upper triangular matrix `U`.

# References 
[1] : Heath algorithm 2.3
"""
function gaussian_elimination!(A)
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
    gaussian_elimination(A)

Return `L` and `U` factorized matrices of `A` explicitly via OOP operation.
"""
function gaussian_elimination(A) 
    A_ = copy(A)
    gaussian_elimination!(A_)
    # Explicitly cast to matrices to show zeros
    L = Matrix(UnitLowerTriangular(A_))
    U = Matrix(UpperTriangular(A_))
    return L, U
end
