# AQUA-Equation of State
[![License: MIT License](https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square)](https://opensource.org/licenses/MIT) ![Release](https://img.shields.io/github/v/release/mnijh/AQUA?include_prereleases)

![Header](AQUA_title.jpg)

A tabulated equation of state for water based on "AQUA: A Collection of H<sub>2</sub>O Equations of State for Planetary Models" by Haldemann et al. (submitted).  It combines the equation of state of Mazevet et al. (2019) with other equation of state at lower temperatures and pressures, namely:

  - Feistel & Wagner (2006) (Ice-Ih)
  - Journaux et al. (2020)  (Ice-II, -III, -V, -VI)
  - French & Redmer (2015)  (Ice VII, X)
  - Wagner & Pruss (2002)   (liquid and gas, below 1 GPa and below 1200 K)
  - Brown (2018)            (liquid and supercritical above 1 GPa, below 6'000 K)
  - CEA package by Gordon (1994) and McBride (1996)

The equation of state was constructed in pressure - temperature space but we also provide density - temperature and density - internal energy tables. Besides pressure and temperature the tables provide data of the density, adiabatic temperature gradient, entropy, internal energy and bulk speed of sound of water. For a limited region where CEA was used, we also provide the mean molecular weight, ionization fraction and dissociation fraction of water.

The repository further contains a simple Fortran module to access the tabulated data. 

This repository is available under the MIT license. If you use content of this repository cite:

```latex
@article{Haldemann_2020a,
  doi = {},
  url = {},
  year  = {submitted},
  month = {},
  volume = {},
  number = {},
  pages = {},
  author = {Haldemann, J., Alibert, Y., Mordasini, Ch. and Benz, W. },
  title = {AQUA: A Collection of $H_2O$ Equations of State for Planetary Models},
  journal = {Astronomy & Astrophysics}
}
```



## Download the Equation of State Tables

The easiest way to get the equation of state tables is cloning this repository. Since the tables are rather large for a `github` repository, we use `git-lfs` (git large file storage) to manage the table files. Hence before cloning the repository please make sure you have installed `git-lfs`, otherwise the tables will show up as empty pointers. Install is pretty easy, 

``` bash
git lfs install
```

will do the job.

After that you can simply clone the repository using

```bash
git clone git@github.com:mnijh/AQUA.git
```

If you don't want to clone the github repository, you can either download the repository as a .zip file or download each table separately, after you opened them here on the github website (see [Tables](Tables)). 

## Using the Fortran module

In the [Fortran](Fortran) directory, we provide a simple Fortran-95 module to evaluate the tabulated data.

In order to use the Fortran module, first clone or download the repository. Next, within your main Fortran program make sure that before you evaluate any AQUA data, to load the tables using one of the load routines:

```Fortran
!! The path variable should point to the directory where the tables are stored.
call LOAD_TABLE_PT(path)
call LOAD_TABLE_RhoT(path)
call LOAD_TABLE_RhoU(path)
```

Then the AQUA tables can be evaluated using one of the functions 

```Fortran
EOS_PT = INTERPOLATE_AQUA_PT(P,T) !! P in Pa, T in K
  
!! EOS_PT contains (density [kg/m^3], adiabatic temperature gradient, entropy [J/kg/K], internal energy [J/kg], bulk speed of sound [m/s], mean molecular weight [g/mol], ionization fraction = N_e/N_tot, dissociation fraction = 1 - N_H2O/N_tot)
```

```Fortran
EOS_RhoT = INTERPOLATE_AQUA_RhoT(Rho,T) !! Rho in kg/m^3, T in K
  
!! EOS_RhoT contains (pressure [Pa], adiabatic temperature gradient, entropy [J/kg/K], internal energy [J/kg], bulk speed of sound [m/s], mean molecular weight [g/mol], ionization fraction = N_e/N_tot, dissociation fraction = 1 - N_H2O/N_tot)
```

```Fortran
EOS_RhoU = INTERPOLATE_AQUA_RhoU(Rho,U) !! !! Rho in kg/m^3, U in J/kg
  
!! EOS_RhoU contains (pressure [Pa], temperature [K], adiabatic temperature gradient, entropy [J/kg/K], bulk speed of sound [m/s], mean molecular weight [g/mol], ionization fraction = N_e/N_tot, dissociation fraction = 1 - N_H2O/N_tot)
```

For the sake of simplicity a bilinear interpolation scheme is used here to evaluate the tabulated data. More accurate interpolation can be achieved using a hermite spline interpolation, e.g. PCHIP - Piecewise Cubic Hermite Interpolant Package.

The following code snippet should illustrate the basic usage as outlined above.

```Fortran
PROGRAM MY_PROGRAM
  USE AQUA_EOS,only: LOAD_TABLE_PT,INTERPOLATE_AQUA_PT,DIM_OUTPUT
  USE AQUA_EOS,only: LOAD_TABLE_RhoT,INTERPOLATE_AQUA_RhoT
  USE AQUA_EOS,only: LOAD_TABLE_RhoU,INTERPOLATE_AQUA_RhoU
  IMPLICIT NONE
  
  CHARACTER(LEN=200) :: path  !! READ_TABLE expects a string of length 200
  REAL(8) :: P,T,EOS_PT(DIM_OUTPUT),EOS_PT(DIM_OUTPUT),EOS_PT(DIM_OUTPUT)
  
  !! The path variable should point to the directory
  !! where the tables are stored.
  path = '../Table'
  !! Before evaluating the first time the equation of state
  CALL LOAD_TABLE_PT(path)
  CALL LOAD_TABLE_RhoT(path)
  CALL LOAD_TABLE_RhoU(path)
  
  P = 1.0D9 !! Pressure in Pa
  T = 999.0 !! Temperature in K
  
  !! Evaluate AQUA P-T Table
  EOS_PT = INTERPOLATE_AQUA_PT(P,T) 
  
  !! EOS_PT contains [density [kg/m^3], adiabatic temperature gradient,
  !! entropy [J/kg/K], internal energy [J/kg], bulk speed of sound [m/s],
  !! mean molecular weight [g/mol], ionization fraction = N_e/N_tot, 
  !! dissociation fraction = 1 - N_H2O/N_tot)
  
  T   = 999.0 !! Temperature in K
  Rho = 1.0D3 !! density [kg/m^3]
  
  !! Evaluate AQUA Rho-T Table
  EOS_RhoT = INTERPOLATE_AQUA_RhoT(Rho,T) 
  
  !! EOS_RhoT contains (pressure [Pa], adiabatic temperature gradient, 
  !! entropy [J/kg/K], internal energy [J/kg], bulk speed of sound [m/s],
  !! mean molecular weight [g/mol], ionization fraction = N_e/N_tot, 
  !! dissociation fraction = 1 - N_H2O/N_tot)
  
  Rho = 1.0D3 !! density [kg/m^3]
  U   = 1.0D6 !! internal energy [J/kg]
  
  !! Evaluate AQUA Rho-U Table
  EOS_RhoU = INTERPOLATE_AQUA_RhoU(Rho,U) 
  
  !! EOS_RhoU contains (pressure [Pa], temperature [K], 
  !! adiabatic temperature gradient, entropy [J/kg/K], bulk speed of sound [m/s], 
  !! mean molecular weight [g/mol], ionization fraction = N_e/N_tot,
  !! dissociation fraction = 1 - N_H2O/N_tot)
  
  RETURN
 END PROGRAM MY_PROGRAM
```
### Compilation
No special dependencies are required to compile the AQUA_EoS module. Using for example the ```gfortan```compiler
```Bash
gfortran -pedantic -c aqua_eos.f95
```
will create the needed files necessary to link to your main code.

### Options
Some options can be set to change the behavior of the Fortran module. The options will mainly activate or deactivate checks if the input is beyond the tabulated range. When set
``` FORTRAN
safe_mode = .true. ![default = .true.]
```
then all options listed below will be evaluated.  If set ``false``, then no checks will be performed. Only deactivate safe mode if you know you will not leave the table domain.
``` FORTRAN
!! If outside tabulated range, write error message
notify_outside_table  = .true. ![default = .true.]
!! If outside tabulated range, stop, otherwise evaluate at last P,T inside tabluated range.
stop_outside_table    = .true. ![default = .true.]
```
The options can be set either in the source code of the module or set them somewhere in the parent program, e.g.,
```FORTRAN
PROGRAM MY_PROGRAM
  USE AQUA_EOS,only: safe_mode,notify_outside_table,stop_outside_table
  IMPLICIT NONE
  
  safe_mode = .true.
  notify_outside_table  = .true.
  stop_outside_table    = .false.
  
  ...
  
END PROGRAM MY_PROGRAM
```

## References

- Haldemann, J., Alibert, Y., Mordasini, C., & Benz, W. submitted

- [Mazevet, S., Licari, A., Chabrier, G., & Potekhin, A. Y. 2019, Astronomy &
  Astrophysics, 621, A128](https://www.aanda.org/articles/aa/full_html/2019/01/aa33963-18/aa33963-18.html)

- [French, M. & Redmer, R. 2015, Physical Review B, 91, 014308](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.91.014308)
- [Wagner, W. & Pruß, A. 2002, Journal of Physical and Chemical Reference Data,
  31, 387](https://aip.scitation.org/doi/10.1063/1.1461829)
- [Feistel, R. & Wagner, W. 2006, Journal of Physical and Chemical Reference
  Data, 35, 1021](https://aip.scitation.org/doi/10.1063/1.2183324)

- [Brown, J. M. 2018, Fluid Phase Equilibria, 463, 18](https://www.sciencedirect.com/science/article/pii/S0378381218300530)
- [Journaux, B., Brown, J. M., Pakhomova, A., et al. 2020, Journal of Geophysical
  Research: Planets, 125, e2019JE006176](https://doi.org/10.1029/2019JE006176)
- McBride, B. J. G. 1996, Computer Program for Calculation of Complex Chem-
  ical Equilibrium Compositions and Applications II. Users Manual and Pro-
  gram Description, Tech. rep., NASA Lewis Research Center
- Gordon, S. 1994, Computer Program for Calculation of Complex Chemical
  Equilibrium Compositions and Applications. Part 1: Analysis, Tech. rep.,
  NASA Lewis Research Center

## License

MIT License

Copyright (c) 2020, Jonas Haldemann

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
