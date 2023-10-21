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

Return the finite difference `mesh`, coefficient matrix `A` and the rhs vector 
`b` given Dirichlet boundary conditions for the left, right, top, and bottom 
of the discretized Laplace equation (i.e., 2D physical domain).

NOTE: The strange indexing arises from the fact that i is a row index for


# References
[1] : Example 11.5 from Heath pg. 461
"""
function laplace_eq_init_arrays(grid_dim_n::Int, lbc, rbc, tbc, bbc)
    # initialize mesh and boundary conditions
    mesh = zeros(grid_dim_n, grid_dim_n)
    mesh[:, 1] .= lbc   # left boundary condition
    mesh[:, end] .= rbc # right boundary condition
    mesh[1, :] .= tbc   # top boundary condition
    mesh[end, :] .= bbc # bottom boundary condition

    # initialize interior point matrices
    n_interior_points_in_one_direction = grid_dim_n - 2
    n_nodes = n_interior_points_in_one_direction^2
    A = zeros(
        n_nodes,
        n_nodes,)
    b = zeros(n_nodes)

    # from approximation at each interior point in mesh from finite difference
    # (Δi, Δj, coefficient)
    stencil = [(0, 0, 4), (1, 0, -1), (-1, 0, -1), (0, 1, -1), (0, -1, -1)]

    # Make the A matrix
    A_row = 1 # flat row index for A matrix
    boundary_node_domain = (0, n_interior_points_in_one_direction+1)

    # for each interior node in the mesh, use the stencil Δi and Δj
    # to determine if u_{i+Δi, j+Δj} is a boundary node... 
    # if it's not, then compute a flat index for the column in which 
    # the coefficient of u_{i+Δi, j+Δj} will be stored
    for j = 1:n_interior_points_in_one_direction 
        for i in 1:n_interior_points_in_one_direction
            for (Δi, Δj, coeff) in stencil
                on_boundary = i+Δi in boundary_node_domain || 
                    j+Δj in boundary_node_domain
                # A is u11, u21, u31 while interior grid is u11, u12, u13
                # which is why flat index here uses j as row index for
                # interior grid and i as column index
                A_col = (j+Δj-1)*n_interior_points_in_one_direction + 
                    i+Δi

                # Update the A and b arrays using the flat index
                if !on_boundary 
                    A[A_row, A_col] = coeff
                else
                    # note: example 11.5 from Heath pg. 461 treats
                    # treats i, j as x, y axticks on mathematical coordinates
                    # left boundary
                    if i+Δi == 0
                        b[A_row] += lbc
                    # right boundary
                    elseif i+Δi == n_interior_points_in_one_direction+1
                        b[A_row] += rbc
                    # bottom boundary
                    elseif j+Δj == 0
                        b[A_row] += bbc
                    # top boundary
                    else
                        b[A_row] += tbc
                    end
                end
            end 
            A_row += 1 # update flat row index
        end
    end 

    return mesh, A, b
end
