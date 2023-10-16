"""
    compute_global_coordinate_array(lx, dims::Int = 1, neles::Int = 1; ly, lz)

Return array of shape [`dims`, `neles+1`]

Need to think about general cases for nn and neles from Simpson2017 since...

# References
[1] : ch. 3.4.1, 5.2, and 6.3 Simpson2017
"""
function compute_global_coordinate_array(
    lx::Int, neles::Int, dims::Int = 1; ly::Int, lz::Int)
    @assert dims in [1, 2, 3] "`dims` must be 1D, 2D, or 3D (i.e., [1, 2, 3])"
    @assert neles > 0 "`neles` must be natural number" 
 
    # initialize global cooordinate/connectivity array
    n_nodes = neles + 1
    g_coord = zeros(Int, dims, n_nodes)

    # compute x-nodal info
    dx = lx ÷ neles # dist between subseq nodes in x-dir
    nx = compute_num_nodes(lx, dx, neles)

    # constants for ind
    x = one_d = 1
    y = two_d = 2
    z = three_d = 3

    # global node counter
    node_counter = 1

    throw("notimplemented")

    # create connectivity arrays based on dimensions of domain
    if dims == one_d
        for i = 1:nx
            g_coord[x, node_counter] = (i - 1)*dx
            node_counter += 1
        end  

    elseif dims == two_d
        dy = 
        ny = 
        for i = 1:nx, j = 1:ny 
            g_coord[x, node_counter] = (i - 1)*dx
            g_coord[y, node_counter] = (j - 1)*dy
            node_counter += 1
        end 

    elseif dims == three_d 
        ny = ly ÷ neles
        nz = lz ÷ neles
        for i = 1:nx, j = 1:ny, z = 1:nz
            
        end        
    end
 
    return g_coord
end

"""
    compute_num_nodes(l_dim, d_dim, neles)

Return number of nodes for a dimension of length `l_dim`, a distance between
nodes in the element `d_dim`, and a number of finite elements `neles`.
Accounts for `neles == 1` since in this case, 2 nodes should be used instead
of just 1 node.
"""
compute_num_nodes(l_dim, d_dim, neles) = (l_dim ÷ d_dim) + (neles == 1 ? 1 : 0) 
