# Multigrid preconditioner.solver for solution to equations arising
# from finite difference discretization of Laplace (?) equation.
# CFD Julia: CFD Julia: A Learning Module Structuring an
# Introductory Course on Computational Fluid Dynamics
# PI Author: https://www.sanlab.org/people

"""
    multigrid!(
        dx, dy, nxf::Int, nyf::Int, nxc::Int, nyc::Int, v1::Int, v2::Int, v3::Int, 
        unf::AbstractMatrix, unc::AbstractMatrix, f::AbstractMatrix,
        ef::AbstractMatrix, ec::AbstractMatrix, r::AbstractMatrix,
        tol = 1e-12)

# Examples 
```julia-repl
julia> # TODO 
```

# Arguments 
- `dx,dy`: grid spacing in x and y direction for fine grid.
- `nxf, nyf`: number of grid points in x and y directions on fine grid.
- `nxc, nyc`: number of grid points in x and y directions on coarse grid
(nxc = nxf/2, nyc = nyf/2).
- `v1, v2, v3`: number of iterations in relaxation step for different grid levels.
- `unf`: numerical solution on fine grid (required solution).
- `unc`: solution to the error correction on coarse grid.
- `f`: source term. 
- `ef`: error correction on fine grid.
- `ec`: error correction on coarse grid. 
- `r`: residual of fine grid. 
- `tol`: Tolerance for convergence. 
"""
function multigrid!(
    dx, dy, nxf::Int, nyf::Int, nxc::Int, nyc::Int, v1::Int, v2::Int, v3::Int, 
    unf::AbstractMatrix, unc::AbstractMatrix, f::AbstractMatrix,
    ef::AbstractMatrix, ec::AbstractMatrix, r::AbstractMatrix,
    tol = 1e-12)
    for k = 1:max_iter

        # call relaxation on fine grid and compute the numerical solution
        gauss_seidel_mg!(nxf, nyf, dx, dy, f, unf, v1)

        # compute the residual and L2 norm
        compute_residual!(nxf, nyf, dx, dy, f, unf, r)
        rms = compute_l2norm(nxf, nyf, r)
        if (rms / init_rms) <= tol 
            break
        end

        # restrict the residual from fine level to coarse level
        restriction!(nxc, nyc, r, ec)

        # set solution zero on coarse grid
        #unc[:, :] = zeros(nxc + 1, nyc + 1)
        unc[:, :] .= 0

        # solve on the coarsest level and relax V3 times
        gauss_seidel_mg!(nxc, nyc, 2.0 * dx, 2.0 * dy, fc, unc, v3)

        # prolongate solution from coarse level to fine level
        prolongation!(nxc, nyc, unc, ef)

        # correct the solution on fine level
        for j = 2:nyf
            for i = 2:nxf
                unf[i, j] = unf[i, j] + ef[i, j]
            end
        end

        # relax v2 times
        gauss_seidel_mg!(nxf, nyf, dx, dy, f, unf, v2)
    end
    return nothing
end

"""
    restriction(nxc, nyc, r, ec)

Project the residuals on the fine grid onto the errors of the coarse grid. 

# Arguments
- `nxf, nyf`: UNUSED. number of grid points in x and y directions on fine grid. 
- `nxc, nyc`: number of grid points in x and y directions on coarse grid
    (nxc = nxf/2, nyc = nyf/2).
- `r`: residual matrix for the Poisson equation on fine grid.
- `ec`: error correction on coarse grid.
"""
function restriction!(nxc, nyc, r, ec)
    for j = 2:nyc
        for i = 2:nxc
            # grid index for fine grid for the same coarse point
            center = 4.0 * r[2*i-1, 2*j-1]
            # E, W, N, S with respect to coarse grid point in fine grid
            grid = 2.0 * (r[2*i-1, 2*j-1+1] + r[2*i-1, 2*j-1-1] +
                          r[2*i-1+1, 2*j-1] + r[2*i-1-1, 2*j-1])
            # NE, NW, SE, SW with respect to coarse grid point in fine grid
            corner = 1.0 * (r[2*i-1+1, 2*j-1+1] + r[2*i-1+1, 2*j-1-1] +
                            r[2*i-1-1, 2*j-1+1] + r[2*i-1-1, 2*j-1-1])
            # restriction using trapezoidal rule
            ec[i, j] = (center + grid + corner) / 16.0
        end
    end
    return nothing
end

"""
    prolongation(nxc, nyc, unc, ef)

Project the errors on the the coarse grid onto the fine grid.

# Arguments 
- `nxf, nyf`: UNUSED. number of grid points in x and y directions on fine grid.
- `nxc, nyc`: number of grid points in x and y directions on coarse grid
(nxc = nxf/2, nyc = nyf/2).
- `unc`: solution to the error correction on coarse grid.
- `ef`: error correction on fine grid.
"""
function prolongation!(nxc, nyc, unc, ef)
    for j = 1:nyc
        for i = 1:nxc
            # direct injection at center point
            ef[2*i-1, 2*j-1] = unc[i, j]
            # east neighnour on fine grid corresponding to coarse grid point
            ef[2*i-1, 2*j-1+1] = 0.5 * (unc[i, j] + unc[i, j+1])
            # north neighbout on fine grid corresponding to coarse grid point
            ef[2*i-1+1, 2*j-1] = 0.5 * (unc[i, j] + unc[i+1, j])
            # north-east neighbour on fine grid corresponding to coarse grid point
            ef[2*i-1+1, 2*j-1+1] = 0.25 * (unc[i, j] + unc[i, j+1] +
                                           unc[i+1, j] + unc[i+1, j+1])
        end
    end
    return nothing
end

"""
    gauss_seidel_mg!(nx, ny, dx, dy, f, un, v)

Gauss-Seidel iterative solver (i.e., smoother, relaxation operator).

# Arguments 
- `nx, ny`: number of grid points in x and y directions.
- `dx, dy`: grid spacing in x and y directions.
- `un`: relaxation solution after fixed number of iterations.
- `f`: source term.
- `v`: number of iterations.
"""
function gauss_seidel_mg!(nx, ny, dx, dy, f, un, v)
    rt = zeros(Float64, nx + 1, ny + 1) # temporary variable
    den = -2.0 / dx^2 - 2.0 / dy^2
    for iteration_count = 1:v
        for j = 2:nx
            for i = 2:ny
                rt[i, j] = f[i, j] - (un[i+1, j] - 2.0 * un[i, j] + un[i-1, j]) / dx^2
                -(un[i, j+1] - 2.0 * un[i, j] + un[i, j-1]) / dy^2
                un[i, j] = un[i, j] + rt[i, j] / den
            end
        end
    end
    return nothing
end

"""
    compute_residual!(
        nx::Int, ny::Int, dx, dy,
        f::AbstractMatrix, un::AbstractMatrix, r::AbstractMatrix)

Compute the residual at k-th step using the new solution
"""
function compute_residual!(
    nx::Int, ny::Int, dx, dy,
    f::AbstractMatrix, un::AbstractMatrix, r::AbstractMatrix)
    for j = 2:ny
        for i = 2:nx
            d2udx2 = (un[i+1, j] - 2 * un[i, j] + un[i-1, j]) / (dx^2)
            d2udy2 = (un[i, j+1] - 2 * un[i, j] + un[i, j-1]) / (dy^2)
            r[i, j] = f[i, j] - d2udx2 - d2udy2
        end
    end
    return nothing
end

"""
    compute_l2norm(mx, ny, r::AbstractMatrix)

Calculate the L-2 norm of the residual vector
"""
function compute_l2norm(nx, ny, r::AbstractMatrix)
    rms = 0.0
    for j = 2:ny
        for i = 2:nx
            rms = rms + r[i, j]^2
        end
    end
    rms = sqrt(rms / ((nx - 1) * (ny - 1)))
    return rms 
end