# SyEM example
Implementation of the divergence-free Synthetic Eddy Method in [Nek5000](https://nek5000.mcs.anl.gov/).

This repository is an example applied to a straight pipe.

Developed by<br />
Jacopo Canton - jcanton@mech.kth.se<br />
Lorenz Hufnagel - hufnagel@kth.se<br />
Philipp Schlatter - pschlatt@mech.kth.se<br />
and based on a code by Oana Marin and Elia Merzari (ANL).

## Documentation
A `setup.py` script is provided to get you started.

The [pipeMeshNek](https://github.com/jcanton/pipeMeshNek) software (automatically fetched by `setup.py`) is used for the generation of the mesh and .rea file parameters.
The svn version 1093 of [Nek5000](https://nek5000.mcs.anl.gov/) (automatically fetched) is used to run the example.

### Setup and execution
 - run `./setup.py` to fetch and compile the necessary tools, including [pipeMeshNek](https://github.com/jcanton/pipeMeshNek) and [Nek5000](https://nek5000.mcs.anl.gov/).
 - to launch the simulation run `mpirun -np X ./nek5000 2>&1 | tee logfile`, where `X` is the number of processes you wish to use (or submit to a queuing system with an appropriate script).

### DF-SyEM
Turbulence at inlet is "generated" with divergence-free synthetic eddy method (DF-SyEMiso)

References:
 - Poletto, Ruggero: Divergence free development of the synthetic eddy method in order to improve synthetic turbulence for embedded les simulations, PhD Thesis, The University of Manchester, 2015.
 - JARRIN, Nicolas, et al. A synthetic-eddy-method for generating inflow conditions for large-eddy simulations. International Journal of Heat and Fluid Flow, 2006, 27. Jg., Nr. 4, S. 585-593.
