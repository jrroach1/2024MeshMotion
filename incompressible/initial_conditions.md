
### Initial conditions for nearly incompressible simulations
*(Updated April 5th, 2023)*

For an ideal gas the sound speed, internal energy, and Mach number are related to the primitive variables via
```math
c = \sqrt{\gamma \frac{p}{\rho}}, \quad e = \frac{1}{\gamma - 1} \frac{p}{\rho}, \quad M = \frac{u}{c}.
```
After some algebra to eliminate the sound speed, 
```math
    e = \frac{u^2}{\gamma(\gamma - 1) M^2}, \quad\text{or}\quad M = \sqrt{\frac{\gamma(\gamma - 1)e}{u^2}}.
```
While $u$ varies throughout the simulation, we assume here that $e$ is roughly constant, so that the maximum Mach number can be obtained from the initial internal energy and the maximum velocity. The latter can be estimated by running an incompressible simulation, which yields maximum velocity magnitudes of $u_1 = 2.10$ and $u_2 = 4.90$ for cylinder cases one and two respectively. For a target Mach number of $M = 0.01$, the resulting initial internal energy values are be
```math
    e_1 = 7.88\times{10}^4 \quad\text{and}\quad e_2 = 4.29\times{10}^5.
```
In this nearly incompressible limit the maximum kinetic energy density of the flow $k = \frac{1}{2} \rho u^2$ is much smaller than the internal energy, so that any change in internal energy due to mechanical or viscous heating from the flow should be negligible. Thus the initial internal energy is a good estimate for the internal energy at any point in the simulation, justifying our earlier assumption.
