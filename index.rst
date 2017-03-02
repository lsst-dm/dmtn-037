



:tocdepth: 1

.. Please do not modify tocdepth; will be fixed when a new Sphinx theme is shipped.

.. note::

   **This technote is not yet published.**

   A description of the algorithm to calculate a model of the true sky and forward model atmospheric effects to create matched templates.

.. sectnum::

Introduction
============

As with any medium, light entering Earth's atmosphere will refract.
If the incident light arrives from any direction other than straight overhead, then refraction will bend its path towards an observer such that it appears to come from a location closer to zenith (see :numref:`fig-atmos_refraction`).
While the bulk refraction is easily corrected during astrometric calibration, water vapor in the atmosphere increases the index of refraction for bluer light, which leads to an angular refraction that is wavelength-dependent.
Since optical telescopes have a finite bandwidth, this means that a point-like source along the direction towards zenith, possibly across multiple pixels.
In :numref:`fig-max_dcr`, I calculate the absolute worst case DCR in each of the LSST bands by calculating the difference in refraction for two laser beams of light from the same position but each tuned to one of the extreme wavelengths of the filter.
The DCR in this scenario is far more severe than a realistic source spectrum, but it illustrates how the effect will scale with wavelength and zenith angle.
Please see `Appendix A: Refraction calculation`_ for details on the calculation of refraction used throughout this note.

.. figure:: /_static/refraction.gif
   :name: fig-atmos_refraction
   :target: http://target.link/url

   Light incident from a source outside the Earth's atmosphere refracts and it's path bends towards zenith.
   To an observer on Earth, the source appears to be at a higher elevation angle :math:`\Delta E` than it truly is.
   For a quick and dirty approximation, :math:`\Delta E` (in arcseconds) :math:`\sim 90 - E` (in degrees) = zenith angle


.. figure:: /_static/DCR_ZA-wavelength.png
   :name: fig-max_dcr
   :target: http://target.link/url

   Calculation of the maximum DCR in each of the LSST filters, as a function of zenith angle. 


DCR iterative forward modeling
==============================

Because refraction is wavelength-dependent :eq:`eqn-refraction`, the smearing of astronomical sources will depend on their intrinsic spectrum. In practice, this appears as an elongation of the measured PSF in an image, and a jitter in the source location leading to mis-subtractions in difference imaging when looking for transient or variable sources.
This effect can be properly corrected when designing a telescope by incorporating an atmospheric dispersion corrector (ADC) in the optical path in front of the detector, but it cannot be removed in software once the photons have been collected.
Instead, we can mitigate the effect by building DCR-matched template images for subtraction, and use our knowledge of DCR to build a deep model of the sky from a collection of images taken under a range of observing conditions. 

In the sections below, I begin by describing the notation and mathematical framework I will use to solve the problem.
I will then lay out the details of my particular implementation of the algorithm, along with a brief discussion of approaches that I have not implemented but might be of interest in the future.
Finally, I will give examples using data simulated using StarFast, as well as images from Decam.

Mathematical framework
----------------------

The point-spread function (PSF) of a point-like source is a combination of many effects arising from the telescope and its environment.
In this note, I will assume that all instrumental effects have been properly measured and accounted for, though this is of course a simplification.
Even without instrumental effects, though, we still have the atmosphere to contend with.
Turbluence in the atmosphere leads to a blurring of images (seeing), and the severity can vary significantly between nights, or even over the same night.
In my initial implementation I will ignore the effect of variable seeing, and use only observations with comparable PSFs, but in order to make full use of all the data available it will have to be addressed in the future (see section #REF).
Thus, for this initial investigation I will assume that the only effect that changes the shape of the PSF over a set of observations is DCR.

:numref:`fig-subband_diagram` below illustrates the approach of this algorithm. Since DCR arises from the change in the index of refraction of the atmosphere across a filter bandpass, if we can build a model sky in smaller sub-bands the effect is greatly reduced.
Further, if DCR is the only effect on the PSF, then we can construct these sub-band models with a small enough bandwidth such that the model to be the same for all observations regardless of airmass and parallactic angle, except for a bulk shift of the entire image.
We only measure the combined image from all sub-bands, however, so those shifted sub-band images must be added together, which results in an apparent elongation of the PSF.

.. figure:: /_static/DCR_subband_diagram.png
   :name: fig-subband_diagram
   :target: http://target.link/url

   A star observed at different parallactic and zenith angles may appear slightly elongated in the zenith direction, but this is due to DCR.
   If the full filter were split into two sub-bands (here **A** and **B**), the star would appear round(er) but shifted slightly from the position measured across the full band.
   If we assume that the 'true' image within each sub-band does not change, then the changing elongation in the **N** observed images along the zenith direction can be attributed to these shifts, which we can calculate with :eq:`eqn-refraction`.
   
Notation
^^^^^^^^

I will use matrix notation throughout this note, with vectors written as lower case letters and 2D matrices written as upper case letters.
In this context, images are written as vectors, with all of the pixel values unwrapped.
To emphasize this point, I have added arrows over the vectors, though these don't convey any additional meaning.
Image data is written as :math:`\overrightarrow{s_i}`, where the subscript :math:`i` loops over the input images, while model data is written as :math:`\overrightarrow{y_\alpha}` with the subscript :math:`\alpha` looping over sub-bands.

The matrix :math:`B_{i\alpha}` encodes the transformation due to DCR of model plane :math:`\alpha` to image :math:`i`, and the reverse transformation is written as :math:`B_{\alpha i}^\star`.
Since the sub-bands have a narrow bandwidth, the effect of DCR is a uniform shift of all pixels, so:

.. math::
   :label: eqn-BB_identity

   B_{\alpha i}^\star B_{i\alpha} = \mathbb{1}

Finally, the measured PSF of each image :math:`i` is given by :math:`Q^{(i)}`, which is a matrix that does not change the size of the image.
Or, to put it in more familiar terms, it represents the convolution of any given image with the measured PSF of image :math:`i`.
Since there is no current motivation to make the PSFs of sub-bands different from each other, one fiducial PSF is used for all models without any index: :math:`P`.

Iterative solution derivation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The image :math:`\overrightarrow{s_i}` is the sum of all of the sub-band models (see :numref:`fig-subband_diagram`), each shifted by the appropriate amount of DCR relative to the effective wavelength of the full filter from :eq:`eqn-DCR`:

.. math::
   :label: eqn-basic_sum

    \sum_\alpha B_{i\alpha}  \overrightarrow{y_\alpha} =  \overrightarrow{s_i}

Applying the reverse shift for one sub-band :math:`\gamma`, we can re-write :eq:`eqn-basic_sum` as:

.. math::
   :label: eqn-iterative_sum

    \overrightarrow{y_\gamma} = B_{\gamma i}^\star  \overrightarrow{s_i} - B_{\gamma i}^\star  \sum_{\alpha  \neq \gamma} B_{i\alpha}  \overrightarrow{y_\alpha}  

While this may not at first appear to help, we can now solve this problem iteratively.
In each iteration, we can solve for a new set of sub-band models :math:`\overrightarrow{y_\gamma}` using the solutions :math:`\overrightarrow{y_\alpha}` from the last iteration as fixed input.

Once we have a set of model :math:`\overrightarrow{y_\gamma}`, we can use that to predict the template for a future observation :math:`k`:

.. math::
   :label: eqn-basic_template

    \parallel \overrightarrow{s_k}\!\!\parallel  = \sum_\alpha B_{k\alpha}  \overrightarrow{y_\alpha}
 

Extension to variable seeing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

While not implemented yet, there is a fairly clear path forward to extend the iterative solution from :eq:`eqn-iterative_sum` to the case where additional effects beyond DCR introduce changes to the PSF
Variable seeing is the primary concern in this case, but in principle instrumental and other effects could be accounted for in this manner as well. 

Now, we need to convolve the model with the measured PSF of the image :math:`Q^{(i)}`, and convolve the image with the fiducial PSF used for the model :math:`P`.
This modifies :eq:`eqn-basic_sum` above:

.. math::
   :label: eqn-psf_sum

   \sum_\alpha B_{i\alpha}  Q^{(i)} \overrightarrow{y_\alpha}  = P  \overrightarrow{s_i} 

Now we can once again apply the reverse shift for one sub-band, and re-write :eq:`eqn-psf_sum` as:

.. math::
   :label: eqn-psf_iterative_sum

    Q^{(i)}\overrightarrow{y_\gamma} = B_{\gamma i}^\star  P \overrightarrow{s_i} - B_{\gamma i}^\star  \sum_{\alpha  \neq \gamma} B_{i\alpha}  Q^{(i)} \overrightarrow{y_\alpha}  

Unfortunately, we now have improved estimates for :math:`Q^{(i)}\overrightarrow{y_\gamma}` when what we really want is :math:`y_\gamma`.
This problem is identical to the standard problem of image co-addition, however, so at this point we would hook in an existing algorithm for combining images with variable PSFs.

Implementation
--------------

There are four main factors to consider when turning :eq:`eqn-iterative_sum` into an effective algorithm: 
what initial solution to use as the starting point for iterations,
what conditioning to apply to the new solution found in each iteration,
how to detect and down-weight contaminated data,
and how to determine when to exit the loop.


Finding the initial solution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Assuming we have no prior spectral information, the best initial guess is that all sub-bands have the same flux in all pixels.
If all model planes are equal at the start, a good guess for the flux distribution within a sub-band is the standard co-add of the input images, divided by the number of model planes being used (since those will be summed).
A proper inverse-variance weighting of the input images as part of the coaddition will help make the best estimate, and if there are many input images we could restrict the coaddition to use only those observed near zenith (with negligible DCR).
An advantage of selecting the simple coadd as the starting point, is that the solution should immediately converge if the input data exhibits no actual DCR effects, such as *i*-band or zenith observations.
However, since this image is only the starting point of an iterative process, the final solution should not be sensitive to small errors at this stage.


Conditioning the iterative solution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A common failure mode of iterative forward modeling algorithms is oscillating solutions.
In these cases, :eq:`eqn-iterative_sum` may produce intermediate solutions for :math:`\overrightarrow{y_\gamma}` with very large amplitude in one iteration, leading to very small amplitude (or negative) solutions in the next iteration, for example.
Conditioning of the solution can mitigate this sort of failure, and also help reach convergence faster.
Some useful types of conditioning include:

* Instead of taking the current solution from :eq:`eqn-iterative_sum` directly, use the average of the current and last solutions.
  This eliminates most instances of oscillating solutions, since it restricts the relative change of the solution between iterations.

* Threshold the solutions.
  Instead of solutions diverging through solutions oscillating between iterations, the solution might 'oscillate' between model planes.
  While it is possible for all of the flux in an image to come from one single model plane, with zero from all others, there are limits.
  Solutions with more flux near a source in a single plane than the initial coadd are likely to be unphysical, and also likely to be paired with deeply negative pixels in the other planes.
  Care must be taken to avoid overly strict thresholds that impair convergence (such as applying the preceding test to even noise-like pixels), but reasonable restrictions can eliminate extreme outliers.

* Frequency regularization.
  In addition to comparing the current solution to the last or initial solutions, we could also apply restrictions on variations between model planes.
  For example, we could calculate the slope (or higher derivatives) of the spectrum for every pixel in the model across the sub-bands, and apply a threshold.
  Any values deviating more than a set amount from the line (or higher order curve) fit by that slope could be fixed to the fit instead, and minimum and maximum slopes could be set.
  While I have written an option within the DCR modeling code to enforce this sort of regularization, in practice I have found the additional benefit to be negligible when combined with the preceding forms of conditioning, and leave it turned off by default.


Weighting the input data
^^^^^^^^^^^^^^^^^^^^^^^^

Weighting of the input data takes two forms:
weights that are static properties of the image (such as the variance plane),
and dynamic weights that may change between iterations.

* Static weights.
  In most cases the static weights will be just the inverse of the variance planes of the images, and best practice is to maintain separate arrays of inverse-variance weighted image values and the corresponding inverse variance values.
  All transformations and convolutions are applied to both equally, and the properly weighted solution is the transformed weighted-image sum divided by the transformed weights sum.

* Dynamic weights.
  The simplest form of dynamic weights is a flag, which indicates whether a particular image is to be used in calculating the new solution with :eq:`eqn-iterative_sum` or not.
  If an estimated template is made for each image using the new solution and equation :eq:`eqn-basic_template`, then those templates should become a better fit to the images with every iteration.
  While it is possible to have a catastrophic failure where the model performs worse for *all* images, it might also improve for most and degrade for a few.
  For example, if there are astrometric errors for one image, the pixel-based model may be misaligned to that image and the subtraction residuals may increase between iterations.
  In that case, that image would hurt the calculation of the overall model more than the additional data was helping it, and that image should be excluded from the next iteration.
  However, in case the apparent divergence was a fluke, convergence should still be tested for that image on all subsequent iterations in case the fit improves with a better model.
  It might be possible to re-calibrate images that are flagged in this way, with the hope that an improved astrometric solution would also improve the fit to the model.

Determining the end condition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Iterative forward modeling does not have an end condition that can be predetermined, and without setting a limit it would run indefinitely.
Possible end conditions include:

* Fixed time / number of iterations.
  The simplest option is to set an upper bound on the number of iterations, to ensure that the loop does exit within a finite time.
  However, the limit should be set high enough that it does not get hit in typical useage.

* Test for convergence.
  There are two types of tests that check for convergence; one that is fast and one that is accurate.
  The fast check simply compares the rate of change of the solution, and if the difference between the new and the last solution is less than a specified threshold fraction of the average of the two the solution has converged.
  The accurate check, on the other hand, creates a template for each image using :eq:`eqn-basic_template`, and calculates a convergence metric from the difference of each image from its template.
  Once the convergence metric changes by less than a specified level between iterations we can safely exit the loop since additional iterations will provide insignificant improvement.
  This is slower since templates must be created for each image, for every iteration, but that has the advantage of allowing convergence to be checked for each individual image at no extra cost (see "Dynamic weights" above).
  If the extra computational cost is not considered prohibitive, then the second test of convergence is far superior, since it is more accurate and enables additional tests and weighting.

* Test for divergence.
  If a template is made for testing convergence, we should naturally test also for divergence.
  If the convergence metric actually degrades with a new iteration, that is a clear sign that some feature of the image is being modeled incorrectly and more iterations will only exacerbate the problem.
  One special case is if only some of the images degrade, because it is possible that they contain astrometric or other errors, and we could choose to continue with those observations flagged (see "Dynamic weights" above).
  Otherwise, it is safest to exit immediately and discard the solution from the current iteration, using the last solution from before it started to diverge instead.

  * One possible modification is to calculate a spatially-varying convergence metric, and mask regions that degrade in future iterations.
    This allows an improved solution to be found even if one very bright feature (such as an improperly-masked saturated star or cosmic ray) is modeled incorrectly.


Note that if a convergence test is used, it should only be allowed to exit the loop after a minimum number of iterations have passed.
It will depend on how the convergence metric is calculated and the choice of initial solution, but the first iteration can show a slight degradation of convergence.


Examples with simulated images
------------------------------


Examples with DECam images
--------------------------


The DCR Sky Model
=================


Simulated source spectra
------------------------

Appendix: Refraction calculation
==================================

While the true density and index of refraction of air varies significantly with altitude, I will follow :cite:`Stone1996` in approximating it as a simple exponential profile in density that depends only on measured surface conditions.
While this is an approximation, it is reportedly accurate to better than 10 milliarcseconds for observations within 65 degrees of zenith, which should be sufficient for normal LSST operations.

The refraction of monochromatic light is given by

.. math::
   :label: eqn-refraction

   R(\lambda) &= r_0 n_0(\lambda) \sin z_0 \int_1^{n_0(\lambda)} \frac{dn}{n \left(r^2n^2 -r_0^2n_0(\lambda)^2\sin^2z_0\right)^{1/2}} \nonumber\\
    &\simeq \kappa (n_0(\lambda) - 1) (1 - \beta) \tan z_0 - \kappa (1 - n_0(\lambda)) \left(\beta - \frac{n_0(\lambda) - 1}{2}\right) \tan^3z_0

where :math:`n_0(\lambda)`, :math:`\kappa`, and :math:`\beta` are given by equations :eq:`eqn-n_lambda`, :eq:`eqn-kappa`, and :eq:`eqn-beta` below. 
The differential refraction relative to a reference wavelength is simply:

.. math::
   :label: eqn-DCR

   \Delta R(\lambda) = R(\lambda) - R(\lambda_{ref})

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
