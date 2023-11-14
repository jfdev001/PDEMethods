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
    laplace_eq_init_arrays(;grid_dim_n=4, lbc=0, rbc=0, tbc=1, bbc=0)

```math
\\begin{equation}
    u_{xx} + u_{yy} = 0,\\ x \\in [0, 1],\\ y \\in [0, 1]
\\end{equation}
```

Return the finite difference `mesh`, coefficient matrix `A` and the rhs vector 
`b` given Dirichlet boundary conditions for the left, right, top, and bottom 
of the discretized Laplace equation (i.e., 2D physical domain).

NOTE: The strange indexing arises from the fact that j is used to index
rows while i is used to index columns in the physical domain (see fig 11.7 [1]),
while j is used to index columns and i is used to index rows in the matrix code.

# Examples
```julia-repl
julia> using PDEMethods: laplace_eq_init_arrays 
julia> mesh, A, b = laplace_eq_init_arrays(grid_dim_n=4, lbc=0, rbc=0, tbc=1, bbc=0)
julia> A \ b
4-element Vector{Float64}:
 0.12499999999999999
 0.12499999999999999
 0.37499999999999994
 0.37499999999999994
```

# References
[1] : Example 11.5 from Heath pg. 461
"""
function laplace_eq_init_arrays(;grid_dim_n::Int=5, lbc=0, rbc=0, tbc=0, bbc=0)
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
                # see note of doc for this justification
                A_col = (j+Δj-1)*n_interior_points_in_one_direction + 
                    i+Δi

                # Update the A and b arrays using the flat index
                if !on_boundary 
                    A[A_row, A_col] = coeff
                else
                    # note: example 11.5 from Heath pg. 461 treats
                    # treats i, j as x, y axticks on mathematical coordinates
                    # @show (A_row, i, j, i+Δi, j+Δj)
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
