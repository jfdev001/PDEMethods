# Domain decomposition techniques using Schwarz methods

using LinearAlgebra: diag


function additive_schwarz()

end 

function restricted_additive_schwarz()

end 

"""
    restriction(N::Int, i::Int, Ω_ids::Vector{Vector{Int}})::Matrix{Bool}

Return boolean matrix used for indexing a domain with `N` dofs
to the `i` domain `Ω` with subdomain indices `Ω_ids`. This function is
fundamentally redundant since one could just as easily do `U[Ω_ids[i]]`; 
however, this function matches writings in papers and textbooks.

# References
[1] : Ch. 1.3. from Dolean2015
"""
function restriction(
    N::Int, i::Int, Ω_ids::Vector{Vector{Int}})::Matrix{Bool}
    Ωᵢ = Ω_ids[i]
    Nᵢ = length(Ωᵢ) 
    Rᵢ = zeros(Bool, Nᵢ, N)
    for (R_row_ix, keep_col_ix) in enumerate(Ωᵢ)
        Rᵢ[R_row_ix, keep_col_ix] = true
    end
    return Rᵢ
end

function restriction(U::Vector, i::Int, Ω_ids::Vector{Vector{Int}}) 
    return restriction(length(N), i, Ω_ids)
end

function extension(U::Vector, i::Int, Ω_ids::Vector{Vector{Int}})
    return transpose(restriction(U, i, Ω_ids))
end

function extension(N::Int, i::Int, Ω_ids::Vector{Vector{Int}})
    return transpose(restriction(N, i, Ω_ids))
end

function extension(Rᵢ::Matrix{Bool})
    return transpose(Rᵢ)
end

function restriction_diagonal(Rᵢ)
    m, n = size(Rᵢ)
    begin_diag_col_ix = findfirst(row -> row == 1, Rᵢ[1, :])
    return Rᵢ[1:m, begin_diag_col_ix:begin_diag_col_ix+m-1]
end  
