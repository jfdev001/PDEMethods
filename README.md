# PDEMethods

A Julia package for implementing numerical methods for canonical PDEs. General notes on PDEs, finite differences, and finite elements are dispersed throughout. Intended for self study only.

## Forward, Centered, Backward Finite Differences

Given some differential equation $u'(x)$, one could approximate this derivative using a Taylor expansion of the function $u(x)$ (soln to diff eq) with derivative $u'(x)$ to some desired order. The general form of the taylor expansion is given below:

$$
f(a) = \sum_{n=0}^N \frac{f^{(n)}(a)}{n!} (x-a)^n
$$

See we are interested in the function perturbed by some infinitesimal $\Delta x$, we approximate
the function $u(x)$ to $N = 2$

$$
u(x + \Delta x) = u(x) + u'(x)(\Delta x) + \frac{1}{2}u''(x)(\Delta x)^2 + \mathcal{O}((\Delta x)^3)
$$

and claim that $x - a = \Delta x$ and since $a$ is just a variable we rename it to $x$ (see [here](https://math.stackexchange.com/questions/254792/how-is-the-taylor-expansion-for-fx-h-derived)). If we truncate the expansion to $N = 1$, then we see that 

$$
u(x + \Delta x) = u(x) + u'(x)(\Delta x) + \mathcal{O}((\Delta x)^2)
$$

and solving for $u'(x)$, we get the **forward difference formula**:

$$
u'(x) \approx \frac{u(x + \Delta x) - u(x)}{\Delta x}
$$

One might also be interested in an infinitesimal perturbation to the "left" of the function at $x$, therefore we perform the same computation for $u(x - \Delta x)$, 

$$
u(x - \Delta x) = u(x) + u'(x)(-\Delta x) + \frac{1}{2}u''(x)(-\Delta x)^2 + \mathcal{O}((-\Delta x)^3)
$$

and obtain the **backward difference formula**

$$
u'(x) \approx \frac{u(x - \Delta x) - u(x)}{-\Delta x} = \frac{u(x) - u(x - \Delta x)}{\Delta x}
$$

We are interested also in the **centered difference formula**, computed below.

$$
\begin{aligned}
u(x + \Delta x) - u(x - \Delta x) &= [u(x) + u'(x)(\Delta x)] - [u(x) + u'(x)(-\Delta x)] \\
&= 2u'(x)(\Delta x) \\
u'(x) \approx \frac{u(x + \Delta x) - u(x - \Delta x)}{2 \Delta x}
\end{aligned}
$$

Since we often need an approximation for the second derivative as well, we add the two series where $N = 2$ together and solve for $u''(x)$ as shown below and get the **central difference formula for the second derivative**:

$$
\begin{aligned}
u(x + \Delta x) + u(x - \Delta x) &= [ u(x) + u'(x)(\Delta x) + \frac{1}{2}u''(x)(\Delta x)^2] \\&+ [u(x) + u'(x)(-\Delta x) + \frac{1}{2}u''(x)(-\Delta x)^2] \\
&= 2u(x) + u''(x)(\Delta x)^2 \\
u''(x) &= \frac{u(x + \Delta x) + u(x - \Delta x) - 2u(x)}{(\Delta x)^2}
\end{aligned} 
$$

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
