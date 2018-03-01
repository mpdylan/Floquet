# Floquet

Floquet theory for discrete Schroedinger operators and orthogonal polynomials on the unit circle

## Floquet theory

Floquet theory is the theory of quasiperiodic solutions to Sturm-Liouville or similar eigenvalue problems with periodic potential. The central object in the theory is the Floquet discriminant, which determines for which values of a spectral parameter (e.g. eigenvalue) such solutions exist.

This package contains routines for calculating Floquet discriminants for three distinct settings:
* Schroedinger operators on the real line
* Discrete Schroedinger operators and Jacobi matrices (related to polynomials orthogonal on the real line)
* Szego recurrence relations (related to polynomials orthogonal on the unit circle)
There is also a submodule OPUC which can compute the relevant orthogonal polynomials, which are used in the

## Requirements
This requires the following Julia packages:
* Polynomials
* ODE
* ApproxFun
* 
