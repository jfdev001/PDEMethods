# Preconditioning algorithms:
# General Form: 
# Requires formulation of `(M^T)Ax = b` where `M` is the preconditioning
# matrix (in this case left side) that decreases the sensitivity of `A`
# to solving for the solutions `x` given perturbation in the load vector
# `b` 
# * Jacobi preconditioner
# * Balancing domain decomposition by constraints (BDDC)
