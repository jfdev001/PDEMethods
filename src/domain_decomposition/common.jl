using BlockDiagonals: BlockDiagonal

"""
    restriction(N::Int, i::Int, Ω_ids::Vector{Vector{Int}})::AbstractMatrix{Bool}

Return boolean matrix used for indexing a domain with `N` dofs
to the `i` domain `Ω` with subdomain indices `Ω_ids`. This function is
fundamentally redundant since one could just as easily do `U[Ω_ids[i]]`; 
however, this function matches writings in papers and textbooks.

# References
[1] : Ch. 1.3. from Dolean2015
"""
function restriction(
    N::Int, i::Int, Ω_ids::Vector{Vector{Int}})::AbstractMatrix{Bool}
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

function extension(Rᵢ::AbstractMatrix{Bool})
    return transpose(Rᵢ)
end

"""
    partition_of_unity_matrices_Ds(Ω_ids::Vector{Vector{Int}}) 

Return a vector of matrices corresponding to the partition of unity.

TODO: How does this generalize to 2D and 3D?

# References
[1] : Ch. 1.3.1 from Dolean2015.
"""
function partition_of_unity_matrices_Ds(Ω_ids::Vector{Vector{Int}})
    # Diagonal matrices such that a 1 indicates "possession" of that dof on Ωᵢ
    Ds = [Matrix(1.0I, length(Ωᵢ), length(Ωᵢ)) for Ωᵢ in Ω_ids]

    # iterate through diagonal matrices and domains 
    # and update the diagonal matrices such that an entry Dᵢⱼ corresponds to
    # 1/(# of times global dof index appears in other subdomains)
    for (D, Ωᵢ) in zip(Ds, Ω_ids)
        for (local_ix, global_ix) in enumerate(Ωᵢ)
            ix_in_n_domains = sum([global_ix in domain for domain in Ω_ids])
            D[local_ix, local_ix] = 1/ix_in_n_domains
        end
    end
    return Ds 
end  

"""
    nicolaide_coarse_space(Ω_ids::Vector{Vector{Int}})::AbstractMatrix

Return `Z` corresponding to the "slow" modes of some preconditioned matrix
for a domain decomposition with partition of unity matrices `Ds`.

# References
[1] : Ch. 4.2 Dolean2015.
"""
function nicolaide_coarse_space(Ω_ids::Vector{Vector{Int}})
    N_dofs = length(Set(vcat(Ω_ids...))) # could be precomputed
    N_subdomains = length(Ω_ids)
    Z  = zeros(N_dofs*2, N_subdomains)
    Ds = partition_of_unity_matrices_Ds(Ω_ids) # could be precomputed
    for i in 1:length(Ω_ids) 
        Rᵢ  = restriction(N_dofs, i, Ω_ids) # could be precomputed
        Nᵢ, N = size(Rᵢ)    # number dofs in Ωᵢ by total number of dofs
        Rᵢᵀ = transpose(Rᵢ) # extension operator
        Dᵢ  = Ds[i]
        Zᵢ  = Rᵢᵀ*Dᵢ*Rᵢ*ones(N)
        Z[i:i+N-1, i] .= Zᵢ
    end
    return Z 
end
