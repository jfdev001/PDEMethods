"""
    compute_global_coordinate_matrix(
        lx, nelex::Int, dims::Int = 1; 
        ly = 0, neley::Int = 0, 
        lz = 0, nelez::Int = 0)

Return matrix of shape [`dims`, `nx` + `ny` + `nz`] where `nx = nelex + 1`,
`ny = neley + 1`, and `nz = nelez + 1` that decribes the global physical
coordinates of the desired domain of dimension `dims`.

# Arguments
- `lx`: Length of physical x-dimension.
- `nelex`: Number of finite elements in x-dimension.
- `dims`: Dimensions of domain.
- `ly = 0`: Same as `lx` but for y-dimension.
- `neley::Int = 0`: Same `nelex` but for y-dimension.
- `lz = 0`: Same as `lx` but for z-dimension.
- `nelez::Int = 0`: Same as `nelex` but for z-dimension

# References
[1] : Ch. 3.4.1, 5.2 (see fig 5.3), and 6.3 of Simpson2017 for conceptual
inspiration, and Ch. 5.6 of Simpson2017 for implementation details.
"""
function compute_global_coordinate_matrix(
    lx, nelex::Int, dims::Int = 1; 
    ly = 0, neley::Int = 0, 
    lz = 0, nelez::Int = 0)
    @assert dims in [1, 2, 3] "`dims` must be 1D, 2D, or 3D (i.e., [1, 2, 3])"

    # initialize global coordinate/connectivity matrix
    g_coord = nothing
 
    # compute number of nodes in each direction
    nx = nelex + 1
    ny = neley + 1
    nz = nelez + 1

    # compute distance between each node in each direction
    dx = lx/nelex
    dy = lx/neley
    dz = lx/nelez

    # constants for indexing and name for dimensions
    x = one_d = 1
    y = two_d = 2
    z = three_d = 3

    # global node counter
    node_counter = 1

    # create connectivity matrix based on dimensions of domain
    if dims == one_d 
        g_coord = zeros(one_d, nx)
        for i = 1:nx
            g_coord[x, node_counter] = (i - 1)*dx
            node_counter += 1
        end  

    elseif dims == two_d
        g_coord = zeros(two_d, nx + ny)
        for i = 1:nx, j = 1:ny 
            g_coord[x, node_counter] = (i - 1)*dx
            g_coord[y, node_counter] = (j - 1)*dy
            node_counter += 1
        end 

    elseif dims == three_d 
        g_coord = zeros(three_d, nx + ny + nz)
        for i = 1:nx, j = 1:ny, k = 1:nz
            g_coord[x, node_counter] = (i - 1)*dx
            g_coord[y, node_counter] = (j - 1)*dy
            g_coord[z, node_counter] = (k - 1)*dz
            node_counter += 1
        end        
    end
 
    return g_coord
end

