# Line Constants

Overhead line constants calculation library in Python.

Functions in the library include:
+ **calc_L_int**: Calculates internal inductance of solid or tubular conductor
+ **calc_GMR**: Calculates geometric mean radius (GMR) of solid or tubular conductor
+ **carsons**: Calculates Carson's earth return correction factors Rp and Xp for self or mutual terms
+ **calc_self_Z**: Calculates self impedance term (in Ohm/km)
+ **calc_mutual_Z**: Calculates mutual impedance term (in Ohm/km)
+ **calc_Dubanton_Z**: Calculates Dubanton approximation for self or mutual impedance (in Ohm/km)  
+ **calc_Z_matrix**: Calculates primitive impedance matrix
+ **calc_Y_matrix**: Calculates primitive admittance matrix
+ **calc_kron_Z**: Calculates Kron reduced matrix

Features currently not supported:
+ Skin effect calculation for internal self impedance
+ Input data validation, e.g. check that all phase conductor vectors are the same size

Validation
----------

Output of line constants library example (simply run line_constants.py in the command line):

![screenshot of example output](/docs/py_output.png?raw=true)

Corresponding line constants calculation in DIgSILENT PowerFactory 2017:

![screenshot of tower geometry in PowerFactory](/docs/PF_geometry.png?raw=true)

Impedance matrix:
![screenshot of PowerFactory impedance matrix](/docs/PF_impedance.png?raw=true)

Admittance matrix:
![screenshot of PowerFactory admittance matrix](/docs/PF_admittance.png?raw=true)

References
----------
+ Dommel, H. W., "Electromagnetic Transients Program Reference Manual (EMTP Theory Book)", Chapter 4, "Overhead Transmission Lines"
+ Arrillaga, J., and Watson, N. R., "Computer Modelling of Electrical Power Systems", 2nd Edition, Wiley, 2005, Chapter 2.6

License & Copyright
===================

Copyright (C) 2017 Julius Susanto. All rights reserved.

The code is distributed under the 3-clause BSD license found in the LICENSE file.
