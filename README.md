# Overview
_PowSpec_EFTofLSS_ performs the calculation of the redshift space power spectra multipoles within the effective field theory of large scale structure.
The serial C++ code performs all the 1-loop computations required for both the real space and redshift space dark matter power spectra. 
A resummation procedure is also implemented to account for the effect of velocity displacements which are poorly controlled in the Standard
Perturbation Theory framework. A useful .ini file format is utilised. That tool is due to *Ben Hoyt* and if you find it useful I encourage you to have a better look [here](https://github.com/benhoyt/inih/). All that is required is the definition of output directories, and specification
of 3 input power spectra, assumed to be in the same format as [CAMB](https://github.com/cmbant/CAMB). The MathematicaSheets/ subdirectory summaries the analytic calculation of the RSD power spectra to one-loop, including the derivation of the required Fabrikant formulae for evaluation of the three-Bessel integrals for the 22-type terms.

Expect to see more bells and whistles as time progresses. Any contributions to the continued development of this package are very welcome. An independent pipeline
implementing the same program by David Seery may be found [here](https://github.com/ds283/LSSEFT).

This version of _PowSpec_EFTofLSS_ has been written by Donough Regan at the University of Sussex.

# Releases
The current release of PowSpec_EFTofLSS is 2017.1. This release can be identified via a DOI linking to a deposit at zenodo.org. That same .tar.gz archives for each release are available directly from GitHub, but for citations please use the zenodo.org versions.

2017.1 (x March 2017) Source code DOI:xxx

# Licensing

PowSpec_EFTofLSS is distributed under the GNU General Public License version 2, or (at your option) any later version. This license is bundled with the source code as LICENSE.txt.

PowSpec_EFTofLSS depends on the GNU GSL Scientific library for the integration and splining routines, and is distributed under the GNU General Public License.
The .ini file reader utilised has been distributed subject to the [new BSD license](https://github.com/benhoyt/inih/blob/master/LICENSE.txt).

# How to cite _PowSpec_EFTofLSS_

Further development of PowSpec_EFTofLSS depends on demonstrating its usefulness to funding agencies. If you use PowSpec_EFTofLSS to produce numerical results supporting your research then we would appreciate a citation to the paper

Renormalization of the matter power spectrum in real and redshift space, Lucía Fonseca de la Bella, Donough Regan and David Seery. arXiv:xxxx.xxxxx DOI:xxx

# Acknowledgments

Development of the PowSpec_EFTofLSS code has been supported by the grant Precision tests of the inflationary scenario, funded by the European Union's Seventh Framework Programme (FP/2007–2013) and ERC Grant Agreement No. 308082.

Some development was additionally funded by the UK Science and Technology Facilities Council via grants ST/I000976/1 and ST/L000652/1, which funded the science programme at the University of Sussex Astronomy Centre from April 2011–March 2014 and April 2014–March 2017, respectively.
