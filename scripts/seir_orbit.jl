using OrdinaryDiffEq
using ProgressMeter
using Plots

# tspan parameters
nyears = 200
days_per_year = 365
total_time = nyears * days_per_year

# Model parameters (see Figure 5.6 description)
mu = 0.02 / days_per_year
sigma = 1/8 # days
gamma = 1/5 # days
omega = 2 * π / days_per_year  # frequency of oscillations per year

# bifurcation parameter
p_values = beta_ones = LinRange(0,0.35,200)

# Initial conditions
beta_zero = 3.6
beta = beta_zero
s0 = (mu  + sigma)*(mu + gamma)/(sigma*beta)
i0 = mu*(1-s0)/(beta*s0)
e0 = (mu + gamma)*i0/sigma
y0 = [s0, e0, i0]

function seir_rule(y, p, t)
    beta_one = p[1]
    s, e, i = y
    beta = beta_zero * (1 + beta_one * cos(omega * t))
  
    sdot = mu - (beta*i + mu)*s
    edot = beta*s*i - (mu  + sigma)*e
    idot = sigma*e - (mu + gamma)*i

    return [sdot, edot, idot]
end

function seir_log_rule(y, p, t)
    y = exp.(y)
    dy = seir_rule(y, p, t)
    return dy ./ y
end     

# compute and plot orbit diagram
b1s = []
inf = []
n_save = 100
@showprogress for beta_one in p_values
    prob = ODEProblem(
        seir_log_rule, log.(y0), (0, total_time), [beta_one])
    X = solve(
      prob, Tsit5(); 
      # save every solution every 365 days (tsteps)
      saveat=days_per_year)
    i = exp.(X[3, :]) # exponentiate infected state var

    # Init conditions for next iteration are 
    # most stable system state (i.e., last vector of solns)
    global y0 = exp.(X[:, end]) 

    # Get the last n-save timesteps since
    # these correspond to steady state solutions
    i_stable = i[length(i)-((n_save)-1):length(i)]

    # update the vectors holding stable infected
    # states and the associated bifurcation parameter
    push!(inf, i_stable...)
    push!(b1s, repeat([beta_one],n_save)...)
end

scatter(
    b1s, inf, 
    markersize=1, alpha=1, 
    color=:black, legend = false,
    yaxis = :log)
