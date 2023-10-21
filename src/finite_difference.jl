# Finite difference schemes
# see Driscoll 5.4 and Heath 8.6
function centered_first_deriv_wrt_pos(U, i, j, h) 
    (U[i+1,j] - U[i-1,j])/(2*h)
end 

function centered_first_deriv_wrt_time(U, i, j, h)
    (U[i,j+1] - U[i,j-1])/(2*h) 
end 

function centered_second_deriv_wrt_pos(U, i, j, h)
    (U[i+1,j] + U[i-1,j] - 2*U[i,j])/h^2
end 

function centered_second_deriv_wrt_time(U, i, j, h) 
    (U[i,j+1] + U[i,j-1] - 2*U[i,j])/h^2
end 

"""
    laplace_eq_init_arrays(grid_dim_n, lbc, rbc, tbc, bbc)

```math
\\begin{equation}
    u_{xx} + u_{yy} = 0,\\ x \\in [0, 1],\\ y \\in [0, 1]
\\end{equation}
```

Return the finite difference mesh, coefficient matrix `A` and the rhs vector 
`b` given Dirichlet boundary conditions for the left, right, top, and bottom 
of the discretized Laplace equation
"""
function laplace_eq_init_arrays(grid_dim_n::Int, lbc, rbc, tbc, bbc)
    # initialize mesh and boundary conditions
    mesh = zeros(grid_dim_n, grid_dim_n)
    mesh[:, 1] .= lbc   # left boundary condition
    mesh[:, end] .= rbc # right boundary condition
    mesh[1, :] .= tbc   # top boundary condition
    mesh[end, :] .= bbc # bottom boundary condition

    # initialize interior point matrices
    n_finite_diff_coeffs = 4
    n_interior_points_in_one_direction = grid_dim_n - 2
    n_nodes = n_interior_points_in_one_direction^2
    A = zeros(
        n_nodes,
        n_finite_diff_coeffs,)
    b = zeros(n_interior_points_in_one_direction)

    # from approximation at each interior point in mesh from finite difference
    stencil = [(0, 0), (-1, 0), (1, 0), (0, 1)]

    # For each interior node of the mesh, build the corresponding
    # row in the A matrix
    for j = 2:grid_dim_n - 1
        for i in 2:grid_dim_n - 1    
            for (A_col_ix, (Δi, Δj)) in enumerate(stencil)
                on_boundary = i+Δi in (1, grid_dim_n) || j+Δj in (1, grid_dim_n)

                # determine coefficient for A matrix
                if !on_boundary
                    if Δi == 0 && Δj == 0
                        A[i-1, A_col_ix] = 4.0
                    else
                        A[i-1, A_col_ix] = -1.0
                    end
                end    
            end
        end
    end    

    return mesh, A, b
end

"""

            for (Δi, Δj, coeff) in stencil
                on_boundary = i+Δi in (1, grid_dim_n) || j+Δj in (1, grid_dim_n)
                if on_boundary 
                    potential_col_boundary = mesh[j, i+Δi]
                    potential_row_boundary = mesh[j+Δj, i]
                    boundary_sum = potential_row_boundary + potential_col_boundary
                    b
                end 
            end 

""" 
