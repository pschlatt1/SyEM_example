# SyEM example
Implementation of the divergence-free Synthetic Eddy Method in [Nek5000](https://nek5000.mcs.anl.gov/).

This repository is an example applied to a straight pipe.

Developed by<br />
Jacopo Canton - jcanton@mech.kth.se<br />
Lorenz Hufnagel - hufnagel@kth.se<br />
Philipp Schlatter - pschlatt@mech.kth.se<br />
and based on a code by Oana Marin and Elia Merzari (ANL).

## Documentation
This branch is ported to v19 of Nek5000. Open issues are (to be worked on) are:
* too many pressure iterations, which means that the velocity field has a non-zero divergence
* restart of the SyEM has been removed for the time being
* parallelisation needs to be looked at
* isotropy is assumed which is not good close to walls.
* statistics are removed.

### DF-SyEM
Turbulence is generated at the inflow boundary with a divergence-free synthetic eddy method (DF-SyEM)

References:
 - Hufnagel, L., Canton, J., Örlü, R., Marin, O., Merzari, E. and Schlatter, P., The three-dimensional structure of swirl-switching in bent pipe flow. Journal of Fluid Mechanics, Vol. 835, pp. 86-101. [doi:10.1017/jfm.2017.749](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/threedimensional-structure-of-swirlswitching-in-bent-pipe-flow/859F4DC8D58983DF53CBEDCDBF7CE8F5#)
 - Hufnagel, L.: On the swirl-switching in developing bent pipe flow with direct numerical simulation [download here](http://kth.diva-portal.org/smash/record.jsf?dswid=4091&pid=diva2%3A963050&c=1&searchType=UNDERGRADUATE&language=en&query=&af=%5B%5D&aq=%5B%5B%7B%22freeText%22%3A%22hufnagel%22%7D%5D%5D&aq2=%5B%5B%5D%5D&aqe=%5B%5D&noOfRows=50&sortOrder=author_sort_asc&sortOrder2=title_sort_asc&onlyFullText=false&sf=all)
 - Poletto, R.: Divergence free development of the synthetic eddy method in order to improve synthetic turbulence for embedded les simulations, PhD Thesis, The University of Manchester, 2015.
 - Jarrin, N., et al. A synthetic-eddy-method for generating inflow conditions for large-eddy simulations. International Journal of Heat and Fluid Flow, 2006, Vol. 27, No. 4, pp. 585-593.
