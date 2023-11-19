module PDEMethods

include("polynomials.jl")
include("quadrature.jl")

include("basis.jl")
include("finite_difference.jl")
include("fem_mesh.jl")

include("direct_methods.jl")
include("iterative_methods.jl")

include("domain_decomposition/common.jl")
include("domain_decomposition/feti.jl")
include("domain_decomposition/bddc.jl")
include("domain_decomposition/schwarz_methods.jl")

end # module PDEMethods
