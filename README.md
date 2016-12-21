# SyEM example
Implementation of the divergence-free Synthetic Eddy Method in [Nek5000](https://nek5000.mcs.anl.gov/).

This repository is an example applied to a straight pipe.

Developed by<br />
Jacopo Canton - jcanton@mech.kth.se<br />
Lorenz Hufnagel - hufnagel@kth.se<br />
Philipp Schlatter - pschlatt@mech.kth.se<br />

## Documentation
A `setup.py` script is provided to get you started.

The [pipeMeshNek](https://github.com/jcanton/pipeMeshNek) software (automatically fetched by `setup.py`) is used for the generation of the mesh and .rea file parameters.
The svn version 1093 of [Nek5000](https://nek5000.mcs.anl.gov/) (automatically fetched) is used to run the example.

### Setup and execution
 - run `./setup.py` to fetch and compile the necessary tools, including [pipeMeshNek](https://github.com/jcanton/pipeMeshNek) and [Nek5000](https://nek5000.mcs.anl.gov/).
 - to launch the simulation run `mpirun -np X ./nek5000 2>&1 | tee logfile`, where `X` is the number of processes you wish to use (or submit to a queuing system with an appropriate script).
