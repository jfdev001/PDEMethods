# Finite difference scheme for solving wave equation
# NOTE: More intelligent backward/forward diff would be given only scalars and 
# then julia broadcasting would be used for vector case 

"""
    forward_euler_first_order_wave(
        scheme::Symbol = :backward; 
        c = 1, xL = 0, xR = 1, nx = 8_000, tf = 2, Δt=0.0001,
        boundary::Symbol = :periodic)

Uses the forward Euler method `y_{k+1} = y_k + Δt*f(t_k, y_k)` for time 
integration along with spatial differencing scheme corresponding to `scheme` as 
either `:backward` or `:forward` to solve the 1D wave equation. The default
parameters lead to an example 1D diffusion equation that whose wave does
has a limited decrease in height over time. Using a coarser grid leads to 
problems here.

TODO: Could stability analysis be used to determine minimum fineness of 
grid needed for stable algorithm? Why are some grids stable and others are
not? 

``
u_t + c u_x = 0 
``

# Examples
```julia-repl
julia> using PDEMethods: forward_euler_first_order_wave
julia> using Plots 
julia> u, x = forward_euler_first_order_wave()
julia> anim = Animation()
julia> for t in 1:size(u,1)÷100:size(u,1)
           plot(x, u[t, :], color=:blue, legend=false, xlims=(x[begin], x[end]))
           frame(anim)
       end 
julia> gif(anim, "assets/waves.gif", fps=15)
```

# References
[1] : ["The First Order Wave Equation" in Solving PDEs](https://aquaulb.github.io/book_solving_pde_mooc/solving_pde_mooc/notebooks/04_PartialDifferentialEquations/04_01_Advection.html)
"""
function forward_euler_first_order_wave(
    scheme::Symbol = :backward; 
    c = 1, xL = 0, xR = 1, nx = 8_000, tf = 2, Δt=0.0001,
    boundary::Symbol = :periodic)

    # these functions only operate on the interior grid points since the 
    # boundary should remain unchanged  
    function backward_diffscheme(u, t) 
        u_t_i = u[t, 2:end-1]       # i=1 to i=n-1, interior only
        u_t_i_min_1 = u[t, 1:end-2] # i=0 to i=n-2, xL and otherwise interior
        r = u_t_i .- u_t_i_min_1 
        return r
    end

    function forward_diffscheme(u, t) 
        u_t_i = u[t, 2:end-1]
        u_t_i_plus_1 = u[t, 3:end]
        r = u_t_i_plus_1 .- u_t_i         
        return r
    end

    # set spatial difference scheme
    diffscheme = if scheme == :backward 
        backward_diffscheme
    elseif scheme == :forward 
        forward_diffscheme
    else
        throw("notimplemented")
    end

    # validate boundary
    @assert boundary ∈ [:dirichlet, :periodic]

    # get rate of change of vars of interst
    Δx = (xR - xL)/(nx-1)

    # get number of timesteps
    nt = Int(tf / Δt) + 1

    # compute domains 
    x = collect(0:Δx:xR)

    # initialize the computational grid
    u = zeros(nt, nx)

    # time initial condition: u(t=0, x)  = ...
    @. u[1, :] = exp(-200*(x-0.25)^2)

    # homogenous dirichlet boundary conditions
    if boundary == :dirichlet
        u[:, 1] .= 0
        u[:, end] .= 0
    end 

    # Forward euler while updating only interior grid points not 
    # subject to boundary conditions above  
    for t in 1:nt-1
        u[t+1, 2:end-1] .= u[t, 2:end-1] .- c*Δt*(diffscheme(u, t))/Δx

        if boundary == :periodic 
            # use backward diff for this bit where i-1 is just end
            u[t+1, end] = u[t, end] - c*Δt*(u[t, end] - u[t, end-1])/Δx
            u[t+1, 1] = u[t+1, end] # periodic right == left 
        end 
    end

    return u, x
end 
