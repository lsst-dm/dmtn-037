



:tocdepth: 1

.. Please do not modify tocdepth; will be fixed when a new Sphinx theme is shipped.

.. note::

   **This technote is not yet published.**

   A description of the algorithm to calculate a model of the true sky and forward model atmospheric effects to create matched templates.

Introduction
============

As with any medium, light entering Earth's atmosphere will refract. If the incident light arrives from any direction other than straight overhead, then refraction will bend its path towards an observer such that it appears to come from a location closer to zenith (see :numref:`fig-atmos_refraction`). While the bulk refraction is easily corrected during astrometric calibration, water vapor in the atmosphere increases the index of refraction for bluer light, which leads to an angular refraction that is wavelength-dependent. Since optical telescopes have a finite bandwidth, this means that a point-like source along the direction towards zenith, possibly across multiple pixels. In :numref:`fig-max_dcr`, I calculate the absolute worst case DCR in each of the LSST bands by calculating the difference in refraction for two laser beams of light from the same position but each tuned to one of the extreme wavelengths of the filter. The DCR in this scenario is far more severe than a realistic source spectrum, but it illustrates how the effect will scale with wavelength and zenith angle. Please see `Appendix A: Refraction calculation`_ for details on the calculation of refraction used throughout this note.

.. figure:: /_static/refraction.gif
   :name: fig-atmos_refraction
   :target: http://target.link/url

   Light incident from a source outside the Earth's atmosphere refracts and it's path bends towards zenith. To an observer on Earth, the source appears to be at a higher elevation angle :math:`\Delta E` than it truly is. For a quick and dirty approximation, :math:`\Delta E` (in arcseconds) :math:`\sim 90 - E` (in degrees) = zenith angle


.. figure:: /_static/DCR_ZA-wavelength.png
   :name: fig-max_dcr
   :target: http://target.link/url

   Calculation of the maximum DCR in each of the LSST filters, as a function of zenith angle. 


DCR Correction
==============

Because refraction is wavelength-dependent :eq:`eqn-refraction`, the smearing of astronomical sources will depend on their intrinsic spectrum. In practice, this appears as an elongation of the measured PSF in an image, and a jitter in the source location leading to mis-subtractions in difference imaging when looking for transient or variable sources. This effect can be properly corrected when designing a telescope by incorporating an atmospheric dispersion corrector (ADC) in the optical path in front of the detector, but it cannot be removed in software once the photons have been collected. Instead, we can mitigate the effect by building DCR-matched template images for subtraction, and use our knowledge of DCR to build a deep model of the sky from a collection of images taken under a range of observing conditions. 

Mathematical framework
----------------------

.. figure:: /_static/DCR_subband_diagram.png
   :name: fig-subband_diagram
   :target: http://target.link/url

   A star observed at different parallactic and zenith angles may appear slightly elongated in the zenith direction, but this is due to DCR. If the full filter were split into two sub-bands (here **A** and **B**), the star would appear round(er) but shifted slightly from the position measured across the full band. If we assume that the 'true' image within each sub-band does not change, then the changing elongation in the **N** observed images along the zenith direction can be attributed to these shifts, which we can calculate with :eq:`eqn-refraction`.
   

.. math::
   :label: eqn-BB_identity

   B_{i\alpha}^\star B_{\alpha i} = \mathbb{1}

.. math::
   :label: eqn-basic_sum

    \sum_\alpha B_{\alpha i}  \overrightarrow{y_\alpha} =  \overrightarrow{s_i}

.. math::
   :label: eqn-iterative_sum

    \overrightarrow{y_\gamma} = B_{i\gamma}^\star  \overrightarrow{s_i} - B_{i\gamma}^\star  \sum_{\alpha  \neq \gamma} B_{\alpha i}  \overrightarrow{y_\alpha}  

.. math::
   :label: eqn-psf_iterative_sum

    Q_i\overrightarrow{y_\gamma} = B_{i\gamma}^\star  P \overrightarrow{s_i} - B_{i\gamma}^\star  \sum_{\alpha  \neq \gamma} B_{\alpha i}  Q_i \overrightarrow{y_\alpha}  

.. math::
   :label: eqn-psf_sum

   \sum_\alpha B_{\alpha i}  Q_i \overrightarrow{y_\alpha}  = P  \overrightarrow{s_i} 

Implementation
--------------


Examples with simulated images
------------------------------


Examples with DECam images
--------------------------


The DCR Sky Model
=================


Simulated source spectra
------------------------

Appendix A: Refraction calculation
==================================

While the true density and index of refraction of air varies significantly with altitude, I will follow :cite:`Stone1996` in approximating it as a simple exponential profile in density, that depends only on measured surface conditions. While this is an approximation, it is reportedly accurate to better than 10 milliarcseconds for observations within 65 degrees of zenith, which should be sufficient for normal LSST operations.

The refraction of monochromatic light is given by

.. math::
   :label: eqn-refraction

   R(\lambda) &= r_0 n_0(\lambda) \sin z_0 \int_1^{n_0(\lambda)} \frac{dn}{n \left(r^2n^2 -r_0^2n_0(\lambda)^2\sin^2z_0\right)^{1/2}} \nonumber\\
    &\simeq \kappa (n_0(\lambda) - 1) (1 - \beta) \tan z_0 - \kappa (1 - n_0(\lambda)) \left(\beta - \frac{n_0(\lambda) - 1}{2}\right) \tan^3z_0

where :math:`n_0(\lambda)`, :math:`\kappa`, and :math:`\beta` are given by equations :eq:`eqn-n_lambda`, :eq:`eqn-kappa`, and :eq:`eqn-beta` below. 

The index of refraction as a function of wavelength :math:`\lambda` (in Angstroms) can be calculated from the relative humidity (:math:`RH`, in percent), surface air temperature (:math:`T`, in Kelvin), and pressure (:math:`P_s` in millibar):

.. math::
   :label: eqn-n_lambda

   n_0( \lambda ) &=\:& 1 + \Delta n_s + \Delta n_w \\
   \\
   \Delta n_s &=\:& \bigg(2371.34 + \frac{683939.7}{130 -\sigma(\lambda)} + \frac{4547.3}{38.9 - \sigma(\lambda)^2}\bigg) D_s \times 10^{-8} \\
   \\
   \Delta n_w &=\:& \big(6487.31 + 58.058 \sigma(\lambda)^2 - 0.71150\sigma(\lambda)^4 + 0.08851\sigma(\lambda)^6\big) D_w \times 10^{-8} \\
   \\
   \sigma(\lambda) &=\:& 10^4/\lambda \;\;\;( \mu m^{-1})
   

Where the density factors for water vapor :math:`D_w` and dry air :math:`D_s` are given by :eq:`eqn-D_w` and :eq:`eqn-D_s` (from  :cite:`Owens67`), and the water vapor pressure :math:`P_w` is calculated from the relative humidity :math:`RH` with :eq:`eqn-P_w`.

.. math::
   :label: eqn-D_w

   D_w = \bigg[1+P_w (1+3.7\times10^{-4}P_w)\bigg(-2.37321\times 10^{-3} + \frac{2.23366}{T} - \frac{710.792}{T^2} + \frac{7.75141\times 10^4}{T^3}\bigg)\bigg] \frac{P_w}{T} 

.. math::
   :label: eqn-D_s

   D_s = \bigg[1 + (P_s - P_w) \bigg( 57.90 \times 10^{-8} -  \frac{9.3250\times 10^{-4}}{T} + \frac{0.25844}{T^2}\bigg)\bigg] \frac{P_s - P_w}{T}

.. math::
   :label: eqn-P_w

   P_w = RH\times 10^{-4}\times e^{(77.3450 + 0.0057 T - 7235.0/T)}/T^{8.2}


The ratio of local gravity at the observing site to :math:`g= 9.81 m/s^2` is given by

.. math::
   :label: eqn-kappa  

    \kappa = g_0/g = 1 + 5.302\times 10^{-3} \sin^2\phi - 5.83\times 10^{-6} \sin^2(2\phi) - 3.15\times 10^{-7} h \label{eqn:kappa}

By assuming an exponential density profile for the atmosphere, the ratio :math:`\beta` of the scale height of the atmosphere to radius of the observing site from the Earth's core can be approximated by:

.. math::
   :label: eqn-beta

   \beta &= \frac{1}{R_\oplus}\int_{0}^\infty \frac{\rho}{\rho_0} dh \nonumber \\
    &\simeq \frac{P_s}{\rho_0g_0 R_\oplus} = \frac{k_BT}{m g_0 R_\oplus} \nonumber \\
    &=  4.5908\times 10^{-6} T 

where :math:`m` is the average mass of molecules in the atmosphere, :math:`R_\oplus` is the radius of the Earth, :math:`k_B` is the Boltzmann constant, and :math:`g_0` is the acceleration due to gravity at the Earth's surface.


References
==========

.. bibliography:: DCR_references.bib
