# MohrCoulomb_StrainSoft
Repo is for incremental driver that is integrated with a strain softening Mohr-Coulomb model. The repo was originally written by A. Yerro-Colom (use [Yerro (2015)](https://upcommons.upc.edu/handle/2117/102412) for reference).

Information about the constitutive relation:
* Uses a smoothed MC yield surface as presented by Abbo and Sloan to alleviate undefined derivates at the kinks
* Has both an Oritiz-Simo integration scheme and an Euler/Newton-Raphson integration scheme
