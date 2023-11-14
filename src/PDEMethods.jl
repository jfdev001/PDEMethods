module PDEMethods

include("polynomials.jl")
include("quadrature.jl")

include("finite_difference.jl")
include("fem_mesh.jl")

include("direct_methods.jl")
include("iterative_methods.jl")

include("domain_decomposition/feti.jl")
include("domain_decomposition/bddc.jl")
include("domain_decomposition/schwarz.jl")

end # module PDEMethods
