# Tabulated EoS to T-ASPP model

This codes fits a tabulated EoS with the T-ASPP model [1].

The tabulated EoS needs to be in cgs units, with columns.

1-rho(g/cm3) 2-G 3-adExp 4-P(g/cm3) 5-eps(ad) 6-cs^2(r) 7-cs^2(c) 8-GruCoeff 

where G is the fundamental derivative, adExp is the adiabatic exponent, eps is the dimensionless internal energy and GruCoeff is the Grüneissen coefficient.

Python codes transforming tabulated EoS from Arizona University database or CompOSE to this format are also included in the repository.

## Libraries
This repository include the following libraries:
* EOS_cold: contains the functions related to the EoS in tabulated or T-ASPP form
* exportFiles: contains the functions to export the model and evaluate the EoS in the high density regime
* findPT: contains the functions used to locate phase transitions in tabulated EoS
* fitPolytropes: contains Levenberg-Marquardt least-squared methods to fit polytropes for given data
* fitModels: contains the fitting polynomials for the phase transitions

## Runing the code

In the main file, one should specify the path of the tabulated EoS, the name of the EoS (the table is assumed in txt format), and the path for the output.

The code is compiled with 

gcc -o taspp main.c EOS_cold.c exportFiles.c findPT.c fitModels.c fitPolytropes.c -lm

Then it can be run. There is an output log through the console specifying the number of phase transitions found and their location, as well as the information of the fitting polytropes and the polynomial models.

## Output of the code.

The code outputs two txt files: eos_TASPP.txt and eos_quant.txt.

The _TASPP file contains all the parameters needed to construct the T-ASPP model. The output is in the format of the eos.par file that the TOV solver and the numeric hydrodynamics in this repository use.

The _quant file contains five columns: 

1-rho  2-G  3-h  4-P  5-csc

These quantities allow for a quick check of the model for being nonconvex, studying the enthalpy, and check the accuracy of the fit of the pressure and classic sound speed.

## Other details

The crust of the EoS is fixed to be the piecewise polytropic fit of SLy. This could be changed by modifying fitPolytropes.c, where it is indicated the kappa and gamma parameters of the last polytrope of the crust.

The fit of the high density polytropes ends at a density of 6e15 g/cm^3. This can be changed to use the whole table at line 282 of fitPolytropes.c

## Phase transitions in the crust

This code is prepared for EoS with phase transitions in the high density region. If there are phase transitions in the crust (rho<2.62780e12) it will probably fail. Minor fixes could be made to the code to ignore phase transitions in the crust and only reconstruct those in the core. Deal with phase transitions in the crust and extend the T-ASPP model there is a work in progress.



[1] M. Berbel and S. Serna, Phys. Rev. D 2023, https://doi.org/10.1103/PhysRevD.108.083031
