# Implements Balancing domain decomposition by constraints to return 
# a preconditioner M to be used in an iterative method (e.g., PCG)
# Zampini2016: "PCBDDC: A CLASS OF DUAL-PRIMAL METHODS IN PETSc"

function bddc()

end

abstract type AbstractDimension end
struct OneD <: AbstractDimension end 
struct TwoD <: AbstractDimension end
struct ThreeD <: AbstractDimension end
struct ND <: AbstractDimension end

"""
    @kwdef mutable struct Dof

Degree of freedom (Dof), i.e., a node in finite element mesh.

TODO: Can a dof satisfy multiple boundary conditions simulataneously?
I.e., a dof satisfies a neumann AND dirichlet boundary condition?
"""
@kwdef mutable struct Dof{dims <: AbstractDimension}
    "Coordinates in global (uniform) mesh"
    global_mesh_coordinates::CartesianIndex 
    "Value of the degree of freedom"
    value::Float64
    "The integer denoting the `i^th` domain(s) to which this dof belongs"
    Ω_is::Vector{Int} 
    "Bool for dof on global boundary ∂Ω"
    on_boundary::Bool 
    "Bool for dof on interface Γ"
    on_interface::Bool 
    "Bool for dof on dirichlet boundary"
    isdirichlet::Bool
    "Bool for dof on neumann boundary"
    isneumann::Bool
    "Bool for dof on robin boundary"
    isrobin::Bool                     
end

@inline function isconnected(x::Dof, y::Dof)::Bool
    norm((x.global_mesh_coordinates - y.global_mesh_coordinates).I) == 1.0
end 

"""
    isedge(x, y, Nx, Ny)::Bool

Return `true` if two different degrees of freedom `x` and `y` with a set `Nx`
of subdomains sharing `x` and a set `Ny` of subdomains sharing `y`, `false`
otherwise.

# References
[1] : Zampini2016 3.1
"""
@inline function isedge(
    x::Dof{TwoD}, y::Dof{TwoD}, Nx::Set{Int}, Ny::Set{Int})::Bool
    length(Nx) == 2 && issetequal(Nx, Ny) && isconnected(x, y) &&
        (!x.isneumann && !x.isdirichlet)
end

"""
    isvertex(x, Nx)::Bool

TODO: How could a dof not be connected to any other dof? Zampini2016 3.1
"""
@inline function isvertex(x::Dof{TwoD}, Nx::Set{Int})::Bool
    length(Nx) == 2 && (!x.isdirichlet || x.isneumann)
end 
