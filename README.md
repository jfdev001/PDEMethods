# PDEMethods

A Julia package for implementing numerical methods for canonical PDEs. General notes on PDEs, finite differences, and finite elements are dispersed throughout. Intended for self study only. The pronoun "we" is used to refer to both the reader and the author (Jared Frazier).

## Forward, Centered, Backward Finite Differences

This section is based on chapter 8 of ref [1].

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

## The Five Point Stencil for the Laplace Equation in 2D

Note that for a PDE like the laplace equation in 2D

$$
\nabla^2f(x, y, t) = 0
$$

if we expand this out, we get

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

and with the Laplace in particular that $\nabla^2f = 0$, we can solve for $f(x, y)$ for a given timestep but given the discretized grid, we claim that the grid is simply a matrix $U$ with elements $u_{i,j}$, so $f(x,y) = u_{i,j}$ and solving appropriately we finally get the following:

$$
u_{i,j} = \frac{u_{i+1,j} + u_{i-1, j} + u_{i, j+1} + u_{i, j-1}}{4}
$$

Since such a problem is an IVP and BVP, we have sufficient information to compute the solutions on the grid using an appropriate iterative method.

## Weak Form of Poisson's Equation in 2D with Dirichlet BCs

This section is based on refs [8] and chapter 3 and 7 of [9].

To derive the weak form of Poisson's equation in 2D, define the following

$$
\begin{equation}
- \nabla \cdot \nabla u = f \text{  Poisson's Equation}  
\end{equation}
$$

with

$$
\begin{equation}
\begin{aligned}
    u(x, 0) &= 0 && \text{Dirichlet Boundaries} \\
    u(1, y) &= 0 \\
    u(x, 1) &= 0 \\
    u(0, y) &= 0 \\
    \text{in } \Omega &\in [0, 1] ^2 && \text{Domain}
\end{aligned}
\end{equation}
$$

Note that in order to derive the weak form, we must multiply equation (1) by $v$ where $v$ is a test function that belongs to a set of functions such that $\forall v \in H_0^1(\Omega)$ where $H_0^1(\Omega)$ is a subset of a Sobolev space $H^1(\Omega)$ elucidated below.

$$
H^1(\Omega) = \{\psi \in C^0(\Omega) \ | \int_{\Omega}(\psi)^2dx < \infty   \}
$$

This reads "$\psi$ is a continuous function on the domain $\Omega$ that is square integrable (that is square integration exists/is finite)".

Then, we are interested in two subsets of $H^1(\Omega)$ that are defined below.

$$
\begin{aligned}
H_E^1(\Omega) &= \{ \psi \in H^1(\Omega)\ |\ \psi\ \text{satisfies all Dirichlet boundary conditions}  \} \\
H_0^1(\Omega) &= \{ \psi \in H^1(\Omega)\ |\ \psi(\vec{x}) = 0\ \text{at all points}\ \vec{x}\ \text{where Dirichlet boundary conditions are specified} \}
\end{aligned}
$$

We thus have a definition for $\forall v \in H_0^1(\Omega)$ that will come up later.

Now, the weak formulation of equation (1) begins with the multiplication by a test function $v$ to both sides of the equation and then integration over the domain $\Omega$.

$$
\begin{equation}
-\int_{\Omega} v \nabla \cdot \nabla u\ d\Omega = \int_{\Omega}vf\ d\Omega
\end{equation}
$$

We want to simplify equation (3) to only have first derivatives (see [here](https://math.stackexchange.com/questions/754511/finite-element-method-weak-formulation?rq=1)) because having second derivatives is more difficult/costly to implement. 

To do this, we define the vector calculus identity

$$
\begin{equation}
\nabla \cdot (v\vec{F}) = v(\nabla \cdot \vec{F}) + \nabla v \cdot \vec{F}
\end{equation} 
$$

where $v$ is a scalar field and $\vec{F}$ is a vector field.

To use this in the derivation of the weak form of Poisson's equation, we simply state that $\vec{F} = \nabla u$, and substituting into equation (4),

$$
\begin{equation}
\nabla \cdot (v \nabla u) = \boxed{v(\nabla \cdot \nabla u)} + \nabla v \cdot \nabla u
\end{equation}
$$

noting that the boxed part of the equation (5) matches the form in equation (3).

Solving equation (5) for $v(\nabla \cdot \nabla u)$, accounting for the negative sign in front of $v$ in equation (3), and then integrating over the domain $\Omega$ gives

$$
\begin{equation}
-\int_{\Omega} v(\nabla \cdot \nabla u)\ d\Omega = -\int_{\Omega} \nabla \cdot (v \nabla u)\ d\Omega + \int_{\Omega} \nabla v \nabla u\ d\Omega 
\end{equation}
$$

and note that the Divergence theorem 

$$
\int_{\Omega} \nabla \cdot \vec{A}\ d\Omega = \oint_{\partial \Omega = \Gamma} \vec{A} \cdot \vec{n}\ d\Gamma
$$

can be applied to 

$$
\begin{equation}
\int_{\Omega} \nabla \cdot (v \nabla u)\ d\Omega = \oint_{\partial \Omega = \Gamma} v \nabla u \cdot \vec{n}\ d\Gamma 
\end{equation}
$$

by assigning $\vec{A} = v \nabla u$. 

Substituting the right hand side equation (7) into equation (6) gives the below equation.

$$
-\int_{\Omega} v(\nabla \cdot \nabla u)\ d\Omega = - \oint_{\partial \Omega = \Gamma} v \nabla u \cdot \vec{n}\ d\Gamma + \int_{\Omega} \nabla v \nabla u\ d\Omega
$$

Finally, since $v(x, y) \in H_0^1(\Omega)$ and thus $v(x, y) = 0$ at all points where Dirichlet boundary conditions are specified, and the boundary domain is defined by $\partial \Omega$, then

$$
\oint_{\partial \Omega = \Gamma} v \nabla u \cdot \vec{n}\ d\Gamma = 0
$$

and the weak form of Poisson's equation is shown below.

$$
\int_{\Omega} \nabla v \nabla u\ d\Omega = \int_{\Omega} v\ d\Omega
$$

# References

[1] Heath, M. T. (2002). Scientific Computing: An Introductory Survey.
Boston: McGraw-Hill. ISBN: 0072399104

[2] Tobin A. Driscoll, Richard J. Braun. Fundamentals of Numerical Computation: Julia Edition. SIAM-Society for
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

[8] Whiteley, J. (2017). Finite Element Methods: A Practical Guide. Springer.

[9] Aerodynamic CFD. "Deriving the Weak Form in 2D for Poisson's Equation".
(2018). url: https://www.youtube.com/watch?v=0dDrzGPFekM

[10] Wikipedia: Vector Calculus Identities. url: https://en.wikipedia.org/wiki/Vector_calculus_identities

[11] Brenner, S.C., Scott, L.R. (2008). "The Mathematical Theory of Finite Element Methods 3ed". Springer.

[12] Wikipedia: Integration by Parts. url: https://en.wikipedia.org/wiki/Integration_by_parts