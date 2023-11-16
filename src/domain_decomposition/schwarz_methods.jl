# Domain decomposition techniques using Schwarz methods

function additive_schwarz()

end 

function restricted_additive_schwarz()

end 

"""
TODO: Doesn't work for all cases...
"""
function restriction(U::Vector, i::Int, Ω_ids::Vector{Vector{Int}})
    N = length(U)
    Ωᵢ = Ω_ids[i]
    Nᵢ = length(Ωᵢ) 
    Rᵢ = zeros(Bool, Nᵢ, N)
    for i in Ωᵢ
        Rᵢ[i, i] = 1
    end
    return Rᵢ
end 

function extension(N::Int, domain_i::Int, domain_indices::Vector{Int})

end
