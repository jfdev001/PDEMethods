
using LinearAlgebra
using SparseArrays

"""
    YourInterpolatorType <: Interpolator{TM <: AbstractMatrix}

Interpolators are composed of a prolongation and restriction linear operator,
i.e., matrix `P` and `R`, respectively, that act to project residuals and
solutions to different levels of the multigrid hierarchy.

# References
[1] : Briggs2000. Chapter 3 and 4 of "A Multigrid Tutorial, 2ed". pp 31-48. 2000.

[2] : AMGCL: Algebraic Multigrid Outline. url: https://amgcl.readthedocs.io/en/latest/amg_overview.html
"""
abstract type Interpolator{TM <: AbstractMatrix} end

struct LinearInterpolator{TM} <: Interpolator{TM} 
    P::TM
    R::TM

    """
        LinearInterpolator(v::Vector)

    Construct LinearInterpolator using the length of the `v` representing the
    number of degrees of freedom at a given level in the multigrid.

    # References
    [1] : [Prolongation matrix formulation](https://stackoverflow.com/questions/78353521/forming-a-matrix-i-e-linear-operator-based-on-a-description-of-coefficients)
    """
    function LinearInterpolator(n_coarse::Int, n_fine::Int = -1)
        if n_fine == -1
            n_fine = 2*n_coarse + 1
        end
        P = zeros(Float64, n_fine, n_coarse)
        P[1, 1] = 1/2

        function veven!(fine_dof, j) 
            if j != 0 && j != n_coarse + 1
                P[fine_dof, j] = 1
            end
        end

        function vodd!(fine_dof, j)
            j != 0 && j != n_coarse+1 && (P[fine_dof, j] = 1/2)
            j+1 != n_coarse+1 && (P[fine_dof, j+1] = 1/2)
        end
   
        j = 1
        n_equations_per_j = 2 
        for fine_dof = 2:n_fine
            fine_dof
            j = fine_dof ÷ n_equations_per_j
            if iseven(fine_dof)
                veven!(fine_dof, j)
            else
                vodd!(fine_dof, j)
            end
        end        
        R = P'
        return new{typeof(P)}(P, R)
    end
end
