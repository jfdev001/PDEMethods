# Direct methods for sparse linear systems 
# should annihilate sparse entries (minimum degree algorithms)
# * LU (Gaussian Elimination)
# * Cholesky Factorization

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
    A_ = A
    for k = 1:n-1
        A_[k, k] == 0 && break
        for i = k+1:n
            multiplier =  A_[i, k]/A_[k, k]
            A_[i, k] = multiplier 
        end 

        for j = k+1:n
            for i = k+1:n
                submatrix_transform = A_[i, j] - A_[i, k]*A_[k, j] 
                A_[i, j] = submatrix_transform 
            end
        end
    end
    return nothing
end 


