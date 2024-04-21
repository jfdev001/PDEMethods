using LinearAlgebra
using SparseArrays

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
julia> A \\ b
4-element Vector{Float64}:
 0.12499999999999999
 0.12499999999999999
 0.37499999999999994
 0.37499999999999994
```

# References
[1] : Example 11.5 from Heath pg. 461
"""
function laplace_eq_init_arrays(;grid_dim_n::Int=4, lbc=0, rbc=0, tbc=1, bbc=0)
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

"""

Return `mesh`, differentiation matrix `A` and rhs `b` (aka f)

```math
\\begin{aligned}
    D_{xy} &= \\frac{-2}{Δx^2} + \\frac{-2}{Δy^2} \\\\
    D_x    &= \\frac{1}{Δx^2} \\\\
    D_y    &= \\frac{1}{Δy^2}
\\end{aligned}
```

Note that the diagonals of the differentiation matrix correspond to the 
entries `u_{ij}` of dofs vector of form

```math
\\begin{bmatrix}
u_{11} \\\\
u_{12} \\\\
\\vdots \\\\
u_{21} \\\\
\\vdots \\\\
u_{N_x, N_y}
\\end{bmatrix},
```

so that differentiation matrix A handles

```math
\\begin{bmatrix}
u_{11}  & u_{12}  & \\cdots \\\\
u_{11}  & u_{12}  & u_{13} & \\cdots \\\\
\\vdots & \\cdots & \\vdots \\\\
\\cdots & \\cdots & u_{N_x, N_y} 
\\end{bmatrix}
```

TODO: Determining Δx and Δy from provided gird Lx, Ly, nx, and ny...
TODO: determine Dx coefficient locations in differentiation matrix. 

# References
[1] : Equation (84) from p. 37 of Pawar2019.

function poisson_eq_init_arrays(
    f = poisson_f; 
    Lxy=1, grid_dim_n=4, lbc=0, rbc=0, tbc=0, bbc=0)
    
    if lbc != 0 || rbc !=0 || tbc != 0 || bbc != 0 
        throw("only 0 bcs are supported currently")
    end

    Δx = Δy = Lxy / grid_dim_n 
    @show Δx

    # initialize mesh and boundary conditions -- btw the mesh is the 
    # same as the solution vector u but includes boundaries 
    mesh = zeros(grid_dim_n, grid_dim_n)
    mesh[:, 1] .= lbc   # left boundary condition
    mesh[:, end] .= rbc # right boundary condition
    mesh[1, :] .= tbc   # top boundary condition
    mesh[end, :] .= bbc # bottom boundary condition
    @show size(mesh)

    # initialize interior point matrices (i.e., differential matrix
    # and rhs vector
    n_interior_points_in_one_direction = grid_dim_n - 2
    n_nodes = n_interior_points_in_one_direction^2
    A = zeros(
        n_nodes,
        n_nodes,)
    b = zeros(n_nodes)
    
    # Construct differential coeffecients from second order central difference
    # of poisson equation 
    D_xy = -2/Δx^2 + -2/Δy^2  # D_xy * u_{ij}
    D_x  = 1/Δx^2             # D_x  * (u_{i+1,j} + u_{i-1, j}) 
    D_y  = 1/Δy^2             # D_y  * (u_{i, j+1} + u_{i, j-1})

    # Get diagonal indices 
    diag_ind = diagind(A, 0)
    sup_diag_ind = diagind(A, 1)
    sub_diag_ind = diagind(A, -1)

    # Set the diagonals u_{ij} and coeffs for u_{j+1,i} and u_{j-1,i}, resp.
    # NOTE: This is the same as sort of formulation as for laplace
    # in the sense that -4uij + 1u{i+1, j} ...
    A[diag_ind] .= D_xy
    A[sup_diag_ind] .= D_y
    A[sub_diag_ind] .= D_y

    @show A

    # Handle u_{i-1, j} and u_{i+1, j} coefficients
    # and populate b
    throw("notimplemented")
    u = mesh 
    A_row = 1 # flat index  
    A_local = [D_y D_xy D_y D_x D_x] 
    for i in 2:n_interior_points_in_one_direction 
        for j in 1:n_interior_points_in_one_direction
            u_ij = u[i, j]
            
            A_row += 1
        end
    end
    
    # Use the mesh to determine how the function `f` evaluates ...
    return mesh, A, b 
end 
"""

"""

Return poisson system of linear equations given a function `f` for the 
right hand side, the domain of x and y, as well as the number of grid points
`nx` (for which `ny = nx = h`). Thus, `grid[nx, ny]` is a corner boundary
(follows ref [2] notation).

# References
[1] : Eq (84) from Pawar2019
[2] : Parts 3 & 4 from https://www.youtube.com/playlist?list=PLcqHTXprNMIN_fw1BIJDNuaX4nX3FHZ9y
"""
function poisson_eq_init_arrays(
    f = poisson_f; 
    xR = 1, xL = 0,
    yR = 1, yL = 0, 
    nx = 5)

    # initialize grid 
    ny = nx
    Δx = (xR - xL) / (nx - 1) 
    Δy = (yR - yL) / (ny - 1)
 
    # wut...
    mesh = Matrix(undef, nx, ny)
    x = xL:Δx:xR
    y = yL:Δy:yR
    for (i, xi) in enumerate(x)
        for (j, yj) in enumerate(y)
            mesh[i, j] = (xi, yj)
        end
    end
    
    # number of interior points for each dimension
    nx_int = nx - 2
    ny_int = ny - 2
    n_eqs = nx_int * ny_int # number of equations/unknowns

    # differentiation matrix 
    D_xy = (-2/Δx^2) + (-2/Δy^2)
    D_x = 1/Δx^2
    D_y = 1/Δy^2

    A = diagm(D_xy*ones(n_eqs)) # u_{ij}
    A[diagind(A, 1)] .= D_y     # u_{i, j+1}
    A[diagind(A, -1)] .= D_y    # u_{i, j-1}

    A[diagind(A, ny - 1)] .= D_x
    A[diagind(A, -(ny - 1))] .= D_x


    return mesh, A 
end 

"""
    poisson_exact_u(x, y)

Return exact solution to poisson equation with [`poisson_f`](@ref) as `f`
evaluated at coordinates `x` and `y`.

# Examples
```julia-repl
julia> using PDEMethods, Plots
julia> x = y = LinRange(0, 1, 512)
julia> u = @. PDEMethods.poisson_exact_u(x', y)
julia> u_mat = zeros(length(x), length(y)) # same result as above
julia> for (i, xi) in enumerate(x)
         for (j, yj) in enumerate(y)
           u_mat[i, j] = PDEMethods.poisson_exact_u(xi, yj)
         end
       end
julia> contourf(x, y, u, levels=15, color=:jet1)
julia> @assert all(u .≈ u_mat)
```

# References
[1] : Equation 85 from Pawar2019.
"""
poisson_exact_u(x, y) = sin(2*π*x)*sin(2*π*y) + 1/16^2*sin(32*π*x)*sin(32*π*y)

"""
    f_on_mesh_grid(f, domain_1D::Vector)::Matrix

Return matrix from function `f` evaluated on a mesh defined by `domain_1D`.

Assumes `f` is of form `f(arg1, arg2)`. 
"""
function f_on_mesh_grid(f, domain_1D::Vector)::Matrix
    n = length(domain_1D)
    M = zeros(n, n)
    for (j, xj) in enumerate(domain_1D)
        for (i, xi) in enumerate(domain_1D)
           M[i, j] = f(xj, xi)
        end
    end
    return M
end 

"""
    poisson_f(x, y)

Return function `f` evaluation for the coordinates `x` and `y`.

# References
[1] : Equation 86 from Pawar2019.
"""
poisson_f(x, y) = -8*π^2*sin(2*π*x)*sin(2*π*y) - 8*π^2*sin(32*π*x)*sin(32*π*y)

## build n x n 2nd central difference matrix
function cdiff2(n)
    h = 2 / (n+1)
    D = spdiagm(-1 => ones(n-1), 0 => -2*ones(n), 1 => ones(n-1))
    D = D / h^2
end

## build n^2 x n^2 discretization of Laplace operator with zero Dirichlet BCs
function lap(n)

    # sparse second difference matrix and identity
    D = cdiff2(n)
    spI = sparse(I,n,n)

    # Kronecker product
    A = kron(spI,D)+kron(D,spI)

end

"""
    init_poisson_problem(n, f = (x, y) -> 5*exp(-10*(x^2 + y^2)))

Return linear system of equations `Ax = b` with initial guess as zeroed `x0`.

# Arguments
- `n`: Number of interior grid points in the x or y direction.
- `f`: Forcing function for poisson equation `-∇u(x,y) = f`
"""
function init_poisson_problem(n, f = (x, y) -> 5*exp(-10*(x^2 + y^2)))
    xfull = collect(LinRange(-1,1,n+2))
    xint = xfull[2:end-1]  
    b = @. f(xint,xint')          # vectorize f evaluated on grid for right-hand side
    b = vec(b)
    A = lap(n)                      # get finite diff discretization of Laplacian
    ndofs = length(b)
    x0 = zeros(ndofs)
    return A, x0, b 
end 

"""
    mit_poisson_problem(n)

Return `A`, `x0`, and `b` for a 2D poisson problem of size `n^2`.
"""
mit_poisson_problem(n) = init_poisson_problem(n)

"""
    poissonSolve(f, n)

Solve Poisson's equation with right-hand side f(x,y) and zero Dirichlet BCs

# Examples
```julia
julia> using PDEMethods: poissonSolve
julia> # solve Poisson's equation with Gaussian right-hand side
julia> n = 100
julia> f(x,y) = 5*exp(-10*(x^2 + y^2)) # could this also be 0??
julia> A, u, b = poissonSolve(f, n)

# References
[1] : https://github.com/mitmath/18303/blob/master/supp_material/poissonFD.ipynb
```
"""
function poissonSolve(f, n)
    xfull = collect(LinRange(-1,1,n+2))
    xint = xfull[2:end-1]  
    b = @. f(xint,xint')          # vectorize f evaluated on grid for right-hand side
    b = vec(b)
    A = lap(n)                      # get finite diff discretization of Laplacian
    @show typeof(A)
    @show typeof(b)
    u = reshape(A \ b, (n,n))       # reshape solution onto n x n grid
    return A, u, b 
end

