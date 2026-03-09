# System Equations

The DGCMG model is based on nonlinear rigid-body dynamics referenced from the attached article.

The function `computeGyroDynamics.m` computes:
- inertia matrix `M(q)`
- Coriolis/centrifugal matrix `C(q,\dot{q})`
- inverse inertia matrix `W = M^{-1}`

The implementation is intended for simulation and control design workflows.
