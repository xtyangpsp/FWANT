# FWANT
Full-wave ambient noise tomography codes and processing scripts, modified from the copy by Yang Shen @ URI and Haiying Gao @ UMASS. 

## Introduction
This package contains codes and scripts that have been heavily modified by Xiaotao Yang based on the original versions with the following major changes:

* Fixed minor bugs in places.
* Added job monitoring scripts for simulation and kernel calculation.
* Added options and functionalities to some processing scripts.
* Seperated some functionalities into standalone functions.
* Improved automated processing.
* Reduced hard-coded values.
* Improved scripts to process large datasets through parallel procedures. These include kernel calculation, kernel assembly, and generation of the inversion matrix. These changes greatly improved the feasibility of handling Tbs of data with several 100s of stations.
* Changed FORTRAN codes for inversion to handle double precission data and 64 bit long integers. The old codes produced wrong results when the number of measurements and the size of kernels exceed a certain number (overflow the integer limit).

## Key references:
<div class="csl-entry">Gao, H., &#38; Shen, Y. (2014). Upper mantle structure of the Cascades from full-wave ambient noise tomography: Evidence for 3D mantle upwelling in the back-arc. <i>Earth and Planetary Science Letters</i>, <i>390</i>, 222–233. https://doi.org/10.1016/j.epsl.2014.01.012</div>

<div class="csl-entry">Shen, Y., Ren, Y., Gao, H., &#38; Savage, B. (2012). An Improved method to extract very-broadband empirical Green’s functions from ambient seismic noise. <i>Bulletin of the Seismological Society of America</i>, <i>102</i>(4), 1872–1877. https://doi.org/10.1785/0120120023</div>

<div class="csl-entry">Yang, X., &#38; Gao, H. (2020). Segmentation of the Aleutian-Alaska Subduction Zone Revealed by Full-Wave Ambient Noise Tomography: Implications for the Along-Strike Variation of Volcanism. <i>Journal of Geophysical Research: Solid Earth</i>, <i>125</i>(11), 1–20. https://doi.org/10.1029/2020JB019677</div>

<div class="csl-entry">Yang, X., &#38; Gao, H. (2018). Full-Wave Seismic Tomography in the Northeastern United States: New Insights Into the Uplift Mechanism of the Adirondack Mountains. <i>Geophysical Research Letters</i>, <i>45</i>(12). https://doi.org/10.1029/2018GL078438</div>

## Feedback
Any questions on this repository or bugs to report, please report throught the "Issues" tab.
