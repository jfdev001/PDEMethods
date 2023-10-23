using FastGaussQuadrature: gausslegendre
using LinearAlgebra: dot

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
    mygausslegendrequad(f, n) 

Gauss-Legendre quadrature of a function `f` for a number of sample points `n`
where `n` also determines the `n^th` Legendre polynomial from which 
quadrature points (roots of `n^th` order Legendre polynomial) `ξₖ` and 
quadrature weights `wₖ` are determined.

# Examples
```julia-repl
julia> using PDEMethods: mygausslegendrequad
julia> f(ξ) = ((1 - ξ)/2)*((1 - ξ)/2)
julia> mygausslegendrequad(f, 2)
0.6666666666666667
```
"""
function mygausslegendrequad(f, n)
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

"""
    gausslegendrequad(f, n)

Return Gauss-Legendre quadrature of function `f` evaluated at `n` integration 
points. See [`mygausslegendrequad`](@ref) for the naive for loop computation.
"""
function gausslegendrequad(f, n)
    legendre_roots_x, weights = gausslegendre(n)
    func_at_legendre_roots_x = f.(legendre_roots_x)
    quadsum = dot(func_at_legendre_roots_x, weights) 
    return quadsum
end 

"""
    gausslegendrequad_lagrange(xs, ys)

This will probably only apply in 2D for quadrilateral elements due to 
presence of xs AND ys (i.e., physical coordinates of mesh).
"""
function gausslegendrequad_lagrange(xs, ys, Lx, Ly)
    # Map xs and ys to expected domain of [-1, 1] X [-1, 1]
    throw("notimplemented")
end 
