# roots of nth legendre polynomial
# https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature
# ch. 4: "Practical Finite Element Modeling in Earth Science using Matlab"
const gauss_quadrature_points::Vector{
    Tuple{Float64, Float64, Vararg{Float64}}} = [  
                                # n-integration points
    (NaN, NaN),                 # 1
    (sqrt(1/3), -sqrt(1/3)),    # 2
    (0., sqrt(3/5), -sqrt(3/5)) # 3
]

const gauss_quadrature_weights::Vector{
    Tuple{Float64, Float64, Vararg{Float64}}} = [ 
                            # n-integration points
    (NaN, NaN),             # 1
    (1., 1.),               # 2
    (8/9, 5/9, 5/9)         # 3 ... should be +/- 5/9 ??
]

"""
    gausslegendrequad(f, n) 

Gauss-Legendre quadrature of a function `f` for a number of sample points `n`
where `n` also determines the `n^th` Legendre polynomial from which 
quadrature points (roots of `n^th` order Legendre polynomial) `ξₖ` and 
quadrature weights `wₖ` are determined.

# Examples
```julia-repl
julia> using PDEMethods: gausslegendrequad
julia> f(ξ) = ((1 - ξ)/2)*((1 - ξ)/2)
julia> gausslegendrequad(f, 2)
0.6666666666666667
```
"""
function gausslegendrequad(f, n)
    @assert n >= 2
    n > 3 && throw("n > 3 not implemented")
    quadsum = 0
    for k = 1:n
        ξₖ = gauss_quadrature_points[n][k]
        wₖ = gauss_quadrature_weights[n][k]
        quadsum += f(ξₖ)*wₖ
    end 
    return quadsum
end 
