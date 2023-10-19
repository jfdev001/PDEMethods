# PDEMethods

A Julia package for implementing numerical methods for canonical PDEs. General notes on PDEs, finite differences, and finite elements are dispersed throughout. Intended for self study only.

## The five point stencil for the laplace equation in 2D

Note that for a PDE like the laplace equation in 2D

$$
\nabla^2f(x, y, t) = 0
$$

if I expand this out, I get

$$
\frac{\partial^2 f}{\partial x^2} + \frac{\partial^2 f}{\partial y^2} = 0
$$

and then using central difference scheme for the second derivative for each term

$$
\begin{aligned}
    \frac{\partial^2 f}{\partial x^2} &= \frac{f(x + \Delta x, y) + f(x - \Delta x, y) - 2f(x, y)}{(\Delta x)^2}\\
    \frac{\partial^2 f}{\partial y^2} &= \frac{f(x, y + \Delta y) + f(x, y - \Delta y) - 2f(x, y)}{(\Delta y)^2}
\end{aligned}
$$

and noting that on 2D grid, i.e. space has been discretized, $\Delta x = \Delta y = h$, then

$$
\begin{aligned}
\frac{\partial^2 f}{\partial x^2} + \frac{\partial^2 f}{\partial y^2} &= \frac{f(x + \Delta x, y) + f(x - \Delta x, y) + f(x, y + \Delta y) + f(x, y - \Delta y) - 4f(x, y)}{h^2} 
\end{aligned}
$$

and with the Laplace in particular that $\nabla^2f = 0$, we can solve for $f(x, y)$ for a given timestep but given the discretized grid, I claim that the grid is simply a matrix $U$ with elements $u_{i,j}$, so $f(x,y) = u_{i,j}$ and solving appropriately I finally ge the following:

$$
u_{i,j} = \frac{u_{i+1,j} + u_{i-1, j} + u_{i, j+1} + u_{i, j-1}}{4}
$$

Since such a problem is an IVP and BVP, we have sufficient information to compute the solutions on the grid using an appropriate iterative method.

# References

[1] Heath, M. T. (2002). Scientific Computing: An Introductory Survey. 
Boston: McGraw-Hill. ISBN: 0072399104 

[2] Chapters 11, 12, and 13 from Tobin A. Driscoll, Richard J. Braun. 
Fundamentals of Numerical Computation: Julia Edition. SIAM-Society for 
Industrial and Applied Mathematics, 2022. url: https://tobydriscoll.net/fnc-julia/frontmatter.html

[3] Georgoulis, M. (2009). Computational Methods for Partial Differential 
Equations. url: http://users.math.uoc.gr/~tsogka/Courses/AESDE-spring2015/Biblio/Georgoulis_notes_new.pdf

[4] Pawar S, San O. CFD Julia: A Learning Module Structuring an Introductory 
Course on Computational Fluid Dynamics. Fluids. 2019; 4(3):159. 
https://doi.org/10.3390/fluids4030159 

[5] Simpson, G. (2017). Practical finite element modeling in Earth science using MATLAB. Wiley. 

[6] Darve, E., Wootters, M (2021). Numerical Linear Algebra with Julia. SIAM.

[7] Recktenwald, G.W. Finite-difference approximations to the heat equation. Mechanical 
Engineering. 2004; 10(1).
