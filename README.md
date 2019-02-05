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
Turbulence is generated at the inflow boundary with a divergence-free synthetic eddy method (DF-SyEM)

References:
 - Hufnagel L., Canton J., Oerlue R., Marin O., Merzari E. and Schlatter P., The three-dimensional structure of swirl-switching in bent pipe flow. Journal of Fluid Mechanics, 835, 86-101. [doi:10.1017/jfm.2017.749](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/threedimensional-structure-of-swirlswitching-in-bent-pipe-flow/859F4DC8D58983DF53CBEDCDBF7CE8F5#)
 - Hufnagel, L.: On the swirl-switching in developing bent pipe flow with direct numerical simulation [download here](http://kth.diva-portal.org/smash/record.jsf?dswid=4091&pid=diva2%3A963050&c=1&searchType=UNDERGRADUATE&language=en&query=&af=%5B%5D&aq=%5B%5B%7B%22freeText%22%3A%22hufnagel%22%7D%5D%5D&aq2=%5B%5B%5D%5D&aqe=%5B%5D&noOfRows=50&sortOrder=author_sort_asc&sortOrder2=title_sort_asc&onlyFullText=false&sf=all)
 - Poletto, R.: Divergence free development of the synthetic eddy method in order to improve synthetic turbulence for embedded les simulations, PhD Thesis, The University of Manchester, 2015.
 - JARRIN, N., et al. A synthetic-eddy-method for generating inflow conditions for large-eddy simulations. International Journal of Heat and Fluid Flow, 2006, 27. Jg., Nr. 4, S. 585-593.
