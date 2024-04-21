module PDEMethods

import AlgebraicMultigrid: _solve

function inspect_amg()
    A, x, b = mit_poisson_problem(10)    
    ml = ruge_stuben(A) 
    soln = _solve(ml, b)
    return A, soln, b
end

export inspect_amg 
export mit_poisson_problem

export Interpolator

include("polynomials.jl")
include("quadrature.jl")

include("basis.jl")
include("finite_difference.jl")
include("fem_mesh.jl")

include("direct_methods.jl")
include("iterative_methods.jl")
include("sparse_smoother.jl")

include("multigrid/interpolator.jl")

include("domain_decomposition/common.jl")
include("domain_decomposition/feti.jl")
include("domain_decomposition/bddc.jl")
include("domain_decomposition/schwarz_methods.jl")

end # module PDEMethods
