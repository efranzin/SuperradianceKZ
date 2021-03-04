[![license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/efranzin/SuperradianceKZ/blob/main/LICENSE)
[![arXiv](https://img.shields.io/badge/arXiv-2102.03152-b31b1b.svg)](https://arxiv.org/abs/2102.03152)

# Amplification factors for bosonic test waves scattered off a Konoplya&ndash;Zhidenko black hole

These data are freely available, but if you make use of them for your research please cite

> [1] Edgardo Franzin, Stefano Liberati, and Mauro Oi, “Superradiance in deformed Kerr black holes” (2021). arXiv: [2102.03152](https://arxiv.org/abs/2102.03152).


## Definitions

For a black hole with mass M, spin parameter a and deformation parameter &eta;,  
the black-hole event horizon r<sub>0</sub> is the largest positive root of r<sup>2</sup> &minus; 2Mr + a<sup>2</sup> &minus; &eta;/r = 0;  
the angular velocity at the horizon is &Omega;<sub>0</sub> = a/(r<sub>0</sub><sup>2</sup> + a<sup>2</sup>) = a/(2Mr<sub>0</sub> + &eta;/r<sub>0</sub>).


## Format of the files

### Massless test fields, s = 0, 1

For a given s, l &ge; s and |m| &le; l.

For each (s, l, m) mode, each `Z_s*_l*_m*.dat` file contains

| black-hole spin | deformation parameter | wave frequency | amplification factor |
|:-:|:-:|:-:|:-:|
| a/M | &eta;/M<sup>3</sup> | &omega; M | Z<sub>slm</sub> |

We give values of a/M = {0.1, &hellip;, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 0.998, 1., 1.01, &hellip;, 1.05, 1.1, 1.15};
for each value of a/M we give values of &eta;/M<sup>3</sup> = {&eta;<sub>&minus;</sub>, &hellip;, 0, 0.01, 0.02, &hellip;, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, &hellip;, 1.}, where &eta;<sub>&minus;</sub> is defined in Eq. (8) of Ref. [1],
and (typically) values of 0 < &omega; < 2|m|&Omega;<sub>0</sub> for m &ne; 0, and 0 < &omega; M < 1 otherwise.

For the non-rotating case, i.e. a = 0, the spectra are degenerate in the azimuthal number m, and for each (s, l) mode, each `Z_a0_s*_l*.dat` file contains

| deformation parameter | wave frequency | absorption factor |
|:-:|:-:|:-:|
| &eta;/M<sup>3</sup> | &omega; M | Z<sub>sl</sub> |

### Massive scalar fields, s = 0, l = m = 1

The file `Z_s0_l1_m1_massive.dat` contains

| black-hole spin | deformation parameter | mass parameter | wave frequency | amplification factor |
|:-:|:-:|:-:|:-:|:-:|
| a/M | &eta;/M<sup>3</sup> | M &mu;<sub>s</sub> | &omega; M | Z<sub>011</sub> |

For now we only give values of a/M = {0.9, 0.95, 0.99} but updates are expected soon.
For each value of a/M we give values of &eta;/M<sup>3</sup> = {&eta;<sub>&minus;</sub>, &hellip;, 0, 0.01, 0.02, &hellip;, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, &hellip;, 1.}, and M &mu;<sub>s</sub> = {0.025, 0.05, &hellip;} as long as &mu;<sub>s</sub> < &Omega;<sub>0</sub>
and (typically) values of &mu;<sub>s</sub> < &omega; < 2 &Omega;<sub>0</sub> &minus; &mu;<sub>s</sub>.


## Data manipulation and plots

The Jupyter notebook [`Plots.ipynb`](Plots.ipynb) provides several plots in addition to those shown in the paper.
It contains plots for the superradiant modes only, for selected values of the black-hole parameters, but it can be easily modified to visualize other modes and/or  other parameter values.


## Notes

These data include the amplification factor for bosonic waves scattered off a Kerr black hole, i.e. with &eta; = 0.


## Acknowledgments

EF acknowledges partial financial support by CNPq Brazil, process no. 301088/2020-9.
EF and SL acknowledge funding from the Italian Ministry of Education and Scientific Research (MIUR) under the grant PRIN MIUR 2017-MB8AEZ.
MO acknowledges partial financial support by the research project “Theoretical and experimental investigations of accreting neutron stars and black holes”, CUP  F71I17000150002, funded by Fondazione di Sardegna.  
The authors thankfully acknowledge Daniele Mura for assistance and computer resources provided by INFN, Sezione di Cagliari.
