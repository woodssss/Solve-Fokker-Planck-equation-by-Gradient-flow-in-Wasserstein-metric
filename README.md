# Fokker Planck equation

In this project, i will introduce a variational scheme for solving Fokker Planck equation, which is known as JKO scheme. Consider a Fokker-Planck equation and we can convert the evolution of the solution into gradient flow in Wasserstein metric. The derivation and explaination are in Explaination.pdf.

## JKO scheme
In fact, the associate gradient flow problem can be solved by JKO scheme.

## Solver
JKO scheme for Fokker-Planck equation is actually a optimization problem with linear constraints. In this project I use Douglas-Rachford algorithm.

# Result and simulation
Given a two bumps initial condition, one can see it eventually converges to global equilibrium.

![grab-landing-page](https://github.com/woodssss/Solve-Fokker-Planck-equation-by-Gradient-flow-in-Wasserstein-metric/blob/master/FP_1D.gif)

## Authors
Wuzhe Xu


