# Numerical integration of differential equations on Lie groups.

This repo provides code to demonstrate the numerical solution of ordinary differential equations (ODEs) defined on Lie groups.  Specialized Lie group integrators are used which preserve the geometric constraint that the Lie-group degrees of freedom are guaranteed to stay on the Lie group manifold (to within machine precision).

The ODE used as a demonstration problem is the so-called *heavy top*, a classic physics problem found in textbooks.  This is a dynamics system where the state variable lies in SO(3) $\times \mathbb{R}^3$: the configuration space is just the space of rotations SO(3), whereas the velocity space is $\mathbb{R}^3$.

Because the repo is intended primarily as a demonstration of the Lie group integrators, for simplicity's sake the code is not generic for any ODE, the way that traditional ODE solvers typically are.  But that is a possible extension.

The Lie Group integrators included in this demonstration are:

**Non-adaptive, constant-timestep integrators**:
* Lie-Euler - The Lie group version of explicit Euler
* Lie-RK2 - The Lie group version of explicit midpoint method, a 2nd-order method
* Lie-RK4 - A variant of Lie-RKMK4 based on reducing the number of commutators
* Lie-RK2CF - A 2nd-order commutator-free method designed for use in a 2(3) adaptive integrator
* Lie-RK3CF - A 3rd-order commutator-free method designed for use in a 2(3) adaptive integrator

**Adaptive-timestep integrators**:
- Lie-RK12 - Adaptive-step method based on the embedded pair of RK2 (midpoint) and explcit Euler
- Lie-RK23CF - Adaptive-step method based on an embedded pair of 2nd and 3rd-order steps, which are *commutator free*.


References
1. A. Iserles, H. Z. Munthe-Kaas, S. P. Norsett, A. Zanna, Lie-group methods, *Acta Numerica* (2005).
2. E. Celledoni, E. Cokaj, A. Leone, D. Murari, B. Owren, Lie group integrators for mechanical systems, *Intl. J. Comp. Math.* (2021).
3. W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery, Numerical Recipes: The Art of Scientific Computing, 3rd Ed. (2007).
