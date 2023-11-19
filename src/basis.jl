# Basis (i.e., shape or test funcs ϕⱼ) functions for FEMs

"""
    ϕ_piecewise_linear(x, j::Int, xs::Vector)

TODO: 2D and 3D.

# Examples
```julia-repl
julia> using PDEMethods, Plots
julia> xs = collect(0:10)
julia> ϕ₁ = [PDEMethods.ϕ_piecewise_linear(x, 1, xs) for x in xs]
julia> ϕ₅ = [PDEMethods.ϕ_piecewise_linear(x, 5, xs) for x in xs]
julia> ϕₖ = [PDEMethods.ϕ_piecewise_linear(x, 11, xs) for x in xs]
julia> p = plot(xs, ϕ₁, label = "ϕ₁")
julia> plot!(p, xs, ϕ₅, label = "ϕⱼ")
julia> plot!(p, xs, ϕₖ, label = "ϕₖ")
```

# References
[1] : Whiteley Ch. 3.3.2.1
"""
function ϕ_piecewise_linear(x, j::Int, xs::Vector)
    k = length(xs) # number of nodes in the mesh
    N = k - 1      # number of elements in mesh
   
    # Internal piecewise linear functions 
    ϕ₁(x, xs) = (xs[1] <= x <= xs[2])*((xs[2] - x)/(xs[2] - xs[1]))
    function ϕⱼ(x, xs) 
        if (xs[j-1] <= x <= xs[j])
            return (x - xs[j-1])/(xs[j] - xs[j-1])
        elseif (xs[j] <= x <= xs[j+1])
            return (xs[j+1] - x)/(xs[j+1] - xs[j])
        end
            return 0
    end
    ϕₖ(x, xs) = (xs[N] <= x <= xs[N+1])*((x - xs[N])/(xs[N+1] - xs[N]))

    # Equations 3.15-3.17 for determining shape func at node
    if (j == 1)
        return ϕ₁(x, xs) 
    elseif (2 <= j <= N)
        return ϕⱼ(x, xs)
    else
        return ϕₖ(x, xs) 
    end
end 
