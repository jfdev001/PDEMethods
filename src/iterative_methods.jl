# Iterative methods for solving linear systems begin with an initial estimate 
# for the solution and successively improve it until the solution is as accurate 
# as desired:
# * jacobi
# * gauss-seidel
# * conjugate gradient method (+ preconditioning?)

"""
    conjugate_gradient_method()

Solve a system of equations `Ax = b` where `A` is a coefficients matrix,
in FEM it would be called the lefthand side matrix consistin of a sum
of a mass matrix and and a stiffness matrix, `x` is a solution vector
to the PDE, and `b` is some constants (i.e., load vector in FEM).
"""
function conjugate_gradient_method()

end 
