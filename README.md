# Rgeogram
An R wrapper to the Geogram library

The Geogram and Graphite libraries are maintained by [Bruno Levy](https://members.loria.fr/blevy/):

http://alice.loria.fr/software/geogram/doc/html/index.html

## Installation

The quickest way to install Rgeogram is via the devtools package:
``` R
install.packages("devtools")
devtools:::install_github("TraME-Project/Rgeogram")
```

Note that Rgeogram requires compilation, so an appropriate development environment is necessary to install the package.
* Windows users should get [Rtools](https://cran.r-project.org/bin/windows/Rtools/). Please ensure that R and Rtools are installed to `C:\` (and not `C:\Program Files`), and that the PATH variables are set correctly (described during the Rtools installation process).
* Mac users should install Xcode and then check [here](https://cran.r-project.org/bin/macosx/tools/) for additional tools (Clang6 and gfortran).

## Example

Solving a 2D problem and recovering the systematic utilities:
``` R
library(Rgeogram)

chi1 <- c(0.2, 0.7, .4, .5)
chi2 <- c(0.1, 0.3, .9, .4)

chi_kj <- cbind(chi1,chi2)/2

res <- otm2D(chi_kj)
otm2D(chi_kj,cbind(rep(0.25,4)))

vtilde <- res$weights

v <- - vtilde + (chi_kj[,1]^2 + chi_kj[,2]^2)
delta_j <- v[4] - v[1:3]
```

## License

* The Rgeogram wrapper code is licensed under the GPL version 2.0. 

* Geogram is licensed under 3-Clause BSD License:

http://alice.loria.fr/software/geogram/doc/html/geogram_license.html

* Graphite is licensed under the GPL:

http://alice.loria.fr/software/graphite/doc/html/graphite_license.html
