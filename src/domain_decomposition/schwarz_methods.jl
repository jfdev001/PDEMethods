# Domain decomposition techniques using Schwarz methods

using LinearAlgebra: diag, I


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

"""
    partition_of_unity_matrices(Ω_ids::Vector{Vector{Int}}) 

Return a vector of matrices corresponding to the partition of unity.

TODO: How does this generalize to 2D and 3D?

# References
[1] : Ch. 1.3.1 from Dolean2015.
"""
function partition_of_unity_matrices(Ω_ids::Vector{Vector{Int}})
    Ds = [Matrix(1.0I, length(Ωᵢ), length(Ωᵢ)) for Ωᵢ in Ω_ids]
    for (D, Ωᵢ) in zip(Ds, Ω_ids)
        for (local_ix, global_ix) in enumerate(Ωᵢ)
            ix_in_n_domains = sum([global_ix in domain for domain in Ω_ids])
            D[local_ix, local_ix] = 1/ix_in_n_domains
        end
    end
    return Ds 
end  
