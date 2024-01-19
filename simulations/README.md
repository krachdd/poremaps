# Benchmarks for Parallel Stokes Solver

### Poiseuille Flow in a Tube 
Compare the velocity in main flow direction $v_3$ to the analytical solution: 
$$v_3(r) = \frac{\Delta p}{4L\mu} (R^2 -r^2)$$ 
See [Batchelor, 1967](https://doi.org/10.1017/CBO9780511800955).

### Rectangular Channel Flow
Compare Flow $Q_i$ to analytical solution, cf. [White and Majdalani, 2006](https://books.google.de/books?id=fl6wPwAACAAJ).

### Sphere-packings 
Compare normalized permeability to Kozeny-Carman equation: 
$$k_1^{KC} = \frac{D^2}{c_{KC}} \frac{\phi^3}{(1-\phi)^2}$$
See [Carman, 1997](https://doi.org/10.1016/S0263-8762(97)80003-2) and [Kozeny, 1927](https://books.google.de/books?id=yERGGwAACAAJ).
